"""
Create protein-ligand system with espaloma.

The system is created using similar schemes found in relative_setup.py from `perses 0.10.1`.
Espaloma is used as the default force field to parameterize both small molecules and proteins
self-consistently.

Note that espaloma and Amber force field is used to first parameterize the small molecule and 
protein, respsectively, to create the solvated complex system. Then the protein from the 
solvated system is extracted and re-parameterized with espaloma. This work around was taken to
handle the issue https://github.com/kntkb/perses/issues/1.

This script was written with following environment: 
- openff-toolkit 0.10.6, openmm 8.0.0, espaloma 0.3.0
"""
import os, sys
import numpy as np
import click
import glob
from rdkit import Chem
import mdtraj as md
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
import openmm.unit as unit
import openmm.app as app
from openmm import MonteCarloBarostat
from openmm import XmlSerializer
#import warnings
import logging
#warnings.filterwarnings("ignore")
logging.basicConfig(filename='logging.log', encoding='utf-8', level=logging.NOTSET)
_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


class CreateSystem(object):
    """
    Create protein-ligand system.

    Notes:
        - Use espaloma force field as default to self-consistently parameterize the protein-ligand system.
        - Use Perses 0.10.1 default parameter settings to setup the system.
    """
    def __init__(self, 
                small_molecule_forcefield = None, 
                forcefield_files = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml'], 
                water_model = 'tip3p', 
                solvent_padding = 9.0 * unit.angstroms, 
                ionic_strength = 0.15 * unit.molar, 
                constraints = app.HBonds, 
                hmass = 3.0 * unit.amu, 
                timestep = 4 * unit.femtoseconds,
                temperature = 300.0 * unit.kelvin, 
                pressure = 1.0 * unit.atmosphere, 
                pme_tol = 2.5e-04, 
                nonbonded_method = app.PME, 
                barostat_period = 50, 
                ):
        self._small_molecule_forcefield = small_molecule_forcefield
        self._forcefield_files = forcefield_files
        self._water_model = water_model
        self._solvent_padding = solvent_padding
        self._ionic_strength = ionic_strength
        self._constraints = constraints
        self._hmass = hmass
        self._timestep = timestep
        self._temperature = temperature
        self._pressure = pressure
        self._pme_tol = pme_tol
        self._nonbonded_method = nonbonded_method
        self._barostat_period = barostat_period


    def create_system(self, protein_file, ligand_file):
        """
        Create protein-ligand system and export serialized system XML file and solvated pdb file.

        Parameters
        ---------
        ligand_file : str
            ligand sdf file. The first ligand entry will be used if multiple ligands are stored.

        Returns
        -------
        None
        """
        _logger.info("Create protein-ligand system")

        # Load ligand
        ext = os.path.splitext(ligand_file)[-1].lower()
        assert ext == '.sdf', f'Ligand file format must be SDF but got {ext}'
        suppl = Chem.SDMolSupplier(ligand_file)
        mols = [x for x in suppl]
        mol = mols[0]
        mol.SetProp("_Name", "MOL")
        offmol = Molecule.from_rdkit(mol)
        ligand_positions = offmol.conformers[0]
        ligand_positions = ligand_positions.in_units_of(unit.nanometers)

        # Load protein pdb
        with open(protein_file, 'r') as f:
            protein = app.PDBFile(f)

        # Join topology
        _logger.info("Merge protein-ligand topology")
        # convert to openmm topology
        protein_topology = protein.topology
        ligand_topology = offmol.to_topology().to_openmm()
        # convert to mdtraj topology
        protein_md_topology = md.Topology.from_openmm(protein_topology)
        ligand_md_topology = md.Topology.from_openmm(ligand_topology)
        # merge topology
        complex_md_topology = protein_md_topology.join(ligand_md_topology)
        complex_topology = complex_md_topology.to_openmm()
        # get number of atoms
        n_atoms_total = complex_topology.getNumAtoms()
        n_atoms_protein = protein_topology.getNumAtoms()
        n_atoms_ligand = ligand_topology.getNumAtoms()
        assert n_atoms_total == n_atoms_protein + n_atoms_ligand, "Number of atoms after merging the protein and ligand topology does not match"
        _logger.info(f"Total atoms: {n_atoms_total} (protein: {n_atoms_protein}, ligand: {n_atoms_ligand}")
        # complex positons
        complex_positions = unit.Quantity(np.zeros([n_atoms_total, 3]), unit=unit.nanometers)
        complex_positions[:n_atoms_protein, :] = protein.positions
        complex_positions[n_atoms_protein:n_atoms_protein+n_atoms_ligand, :] = ligand_positions

        # Initialize system generator
        forcefield_kwargs = {'removeCMMotion': True, 'ewaldErrorTolerance': self._pme_tol, 'constraints' : self._constraints, 'rigidWater': True, 'hydrogenMass' : self._hmass}
        periodic_forcefield_kwargs = {'nonbondedMethod': self._nonbonded_method}
        barostat = MonteCarloBarostat(self._pressure, self._temperature, self._barostat_period)
        self._system_generator = SystemGenerator(
            forcefields=self._forcefield_files, forcefield_kwargs=forcefield_kwargs, periodic_forcefield_kwargs = periodic_forcefield_kwargs, barostat=barostat, 
            small_molecule_forcefield=self._small_molecule_forcefield, molecules=offmol, cache=None)

        # Solvate system
        _logger.info("Solvate system")
        modeller = app.Modeller(complex_topology, complex_positions)
        modeller.addSolvent(self._system_generator.forcefield, model=self._water_model, padding=self._solvent_padding, ionicStrength=self._ionic_strength)
        
        # Get topology and position
        self._solvated_topology = modeller.getTopology()
        self._solvated_positions = modeller.getPositions()
        
        # Create system
        solvated_system = self._system_generator.create_system(self._solvated_topology)
        
        # Save system and export solvated system
        if "espaloma" in self._small_molecule_forcefield:
            app.PDBFile.writeFile(self._solvated_topology, self._solvated_positions, file=open('complex-solvated.pdb', 'w'))
            self._regenerate_espaloma_system()
        else:
            app.PDBFile.writeFile(self._solvated_topology, self._solvated_positions, file=open('complex-solvated.pdb', 'w'))
            # Minimize
            self._minimize(self._solvated_topology, self._solvated_positions, solvated_system)


    def _regenerate_espaloma_system(self):
        """
        Regenerate system with espaloma. Parameterization of protein and ligand self-consistently.

        Reference
        ---------
        https://github.com/kntkb/perses/blob/support-protein-espaloma/perses/app/relative_setup.py#L883

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        _logger.info("Regenerate system with espaloma")

        # Check protein chains
        mdtop = md.Topology.from_openmm(self._solvated_topology)
        chain_indices = [ chain.index for chain in self._solvated_topology.chains() ]
        protein_chain_indices = [ chain_index for chain_index in chain_indices if mdtop.select(f"protein and chainid == {chain_index}").any() ]
        _logger.info(f"Protein chain indices: {protein_chain_indices}")

        # Check conflicting residue names
        conflict_resnames = [ residue.name for residue in mdtop.residues if residue.name.startswith("XX") ]
        if conflict_resnames:
            raise Exception('Found conflict residue name in protein.')

        # Initialize
        self._new_solvated_topology = app.Topology()
        self._new_solvated_topology.setPeriodicBoxVectors(self._solvated_topology.getPeriodicBoxVectors())
        new_atoms = {}

        # Regenerate protein topology
        chain_counter = 0
        _logger.info(f"Regenerating protein topology...")
        for chain in self._solvated_topology.chains():
            new_chain = self._new_solvated_topology.addChain(chain.id)
            # Convert protein into a single residue
            if chain.index in protein_chain_indices:
                resname = f'XX{chain_counter:01d}'
                resid = '1'
                chain_counter += 1
                new_residue = self._new_solvated_topology.addResidue(resname, new_chain, resid)
            for i, residue in enumerate(chain.residues()):
                if residue.chain.index not in protein_chain_indices:
                    new_residue = self._new_solvated_topology.addResidue(residue.name, new_chain, residue.id)
                for atom in residue.atoms():
                    new_atom = self._new_solvated_topology.addAtom(atom.name, atom.element, new_residue, atom.id)
                    new_atoms[atom] = new_atom

        # Regenerate bond information
        for bond in self._solvated_topology.bonds():
            if bond[0] in new_atoms and bond[1] in new_atoms:
                self._new_solvated_topology.addBond(new_atoms[bond[0]], new_atoms[bond[1]])
        
        # Save the updated complex model as pdb
        complex_espaloma_filename = f"complex-solvated-espaloma.pdb"
        if not os.path.exists(complex_espaloma_filename):
            with open(complex_espaloma_filename, 'w') as outfile:
                app.PDBFile.writeFile(self._new_solvated_topology, self._solvated_positions, outfile)
        
        # Seperate proteins into indivdual pdb files according to chain ID.
        protein_espaloma_filenames = glob.glob("protein-espaloma-*.pdb")
        if not protein_espaloma_filenames:
            for chain_index in protein_chain_indices:
                t = md.load_pdb(complex_espaloma_filename)
                indices = t.topology.select(f"chainid == {chain_index}")
                t.atom_slice(indices).save_pdb(f"protein-espaloma-{chain_index}.pdb")
            protein_espaloma_filenames = glob.glob("protein-espaloma-*.pdb")
        
        # Load individual protein structure into openff.toolkit.topology.Molecule
        protein_molecules = [ Molecule.from_file(protein_filename) for protein_filename in protein_espaloma_filenames ]
        
        # We already added small molecules to template generator when we first created ``self._system_generator``.
        # So we only need to add protein molecule to template generator (EspalomaTemplateGenerator).
        self._system_generator.template_generator.add_molecules(protein_molecules)
        # Regenerate system with system generator.
        self._new_solvated_system = self._system_generator.create_system(self._new_solvated_topology)
        # Minimize (use old topology to renumber the atoms and residues)
        #self._minimize(self._new_solvated_topology, self._solvated_positions, self._new_solvated_system)
        self._minimize(self._solvated_topology, self._solvated_positions, self._new_solvated_system)


    def _minimize(self, topology, positions, system):
        """
        Minimize solvated system
        """
        _logger.info(f"Minimize solvated system")

        from openmm import LangevinMiddleIntegrator
        integrator = LangevinMiddleIntegrator(self._temperature, 1/unit.picosecond, self._timestep)
        self._simulation = app.Simulation(topology, system, integrator)
        self._simulation.context.setPositions(positions)
        self._simulation.minimizeEnergy(maxIterations=100)
        self._export_xml(system)


    def _export_xml(self, system):
        _logger.info(f"Serialize and export system")
        state = self._simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
        # Save system
        with open("system.xml", "w") as wf:
            xml = XmlSerializer.serialize(system)
            wf.write(xml)
        # Save and serialize the final state
        with open("state.xml", "w") as wf:
            xml = XmlSerializer.serialize(state)
            wf.write(xml)
        # Save the final state as a PDB
        with open("state.pdb", "w") as wf:
            app.PDBFile.writeFile(
                self._simulation.topology,
                self._simulation.context.getState(
                    getPositions=True,
                    enforcePeriodicBox=False).getPositions(),
                    #enforcePeriodicBox=True).getPositions(),
                    file=wf,
                    keepIds=True
            )
        # Save and serialize integrator
        with open("integrator.xml", "w") as wf:
            xml = XmlSerializer.serialize(self._simulation.integrator)
            wf.write(xml)


@click.command()
@click.option('-p', '--protein_file', required=True, default='protein.pdb', help='protien PDB file')
@click.option('-l', '--ligand_file', required=True, default='ligands.sdf', help='ligand SDF file')
# espaloma-0.3.0.pt is indentical to espaloma-0.3.0rc6.pt and espaloma-0.3.1.pt
@click.option('-f', '--small_molecule_forcefield', default='espaloma-0.3.0.pt', help='Small molecule force field')  
def cli(**kwargs):
    protein_file = kwargs['protein_file']
    ligand_file = kwargs['ligand_file']
    small_molecule_forcefield = kwargs['small_molecule_forcefield']
    
    if 'openff' in small_molecule_forcefield:
        ff_base_path = "/home/takabak/.offxml"
    else:
        ff_base_path = "/home/takabak/.espaloma"

    small_molecule_forcefield = os.path.join(ff_base_path, small_molecule_forcefield)
    if not os.path.exists(small_molecule_forcefield):
        raise FileNotFoundError(f'{small_molecule_forcefield} not found')

    c = CreateSystem(small_molecule_forcefield=small_molecule_forcefield)
    c.create_system(protein_file, ligand_file)


if __name__ == "__main__":
    cli()