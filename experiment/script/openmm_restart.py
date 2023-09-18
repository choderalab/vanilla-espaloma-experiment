import os, sys
import numpy as np
import click
import mdtraj as md
from mdtraj.reporters import NetCDFReporter
from openmm.app import PDBFile, Simulation, CheckpointReporter, StateDataReporter
from openmm import XmlSerializer
#import warnings
#warnings.filterwarnings("ignore")


def export_xml(simulation, system):
    """
    Export system and state
    """
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
    # Save and serialize system
    system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    with open("system.xml", "w") as wf:
        xml = XmlSerializer.serialize(system)
        wf.write(xml)
    # Save and serialize the final state
    with open("state.xml", "w") as wf:
        xml = XmlSerializer.serialize(state)
        wf.write(xml)
    # Save the final state as a PDB
    with open("state.pdb", "w") as wf:
        PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(
                getPositions=True,
                enforcePeriodicBox=False).getPositions(),
                #enforcePeriodicBox=True).getPositions(),
                file=wf,
                keepIds=True
        )
    # Save and serialize integrator
    with open("integrator.xml", "w") as wf:
        xml = XmlSerializer.serialize(simulation.integrator)
        wf.write(xml)


def run(restart_prefix, initial_prefix):
    checkpoint_frequency = 25000  # 10ns
    logging_frequency = 250000  # 1ns
    netcdf_frequency = 250000  # 1ns
    nsteps = 25000000  # 100ns
    
    # Deserialize system file and load system
    with open(os.path.join(restart_prefix, 'system.xml'), 'r') as f:
        system = XmlSerializer.deserialize(f.read())

    # Deserialize integrator file and load integrator
    with open(os.path.join(restart_prefix, 'integrator.xml'), 'r') as f:
        integrator = XmlSerializer.deserialize(f.read())

    # Set up simulation
    pdbfile = os.path.join(initial_prefix, 'complex-solvated.pdb')
    initial_complex = PDBFile(pdbfile)
    simulation = Simulation(initial_complex.topology, system, integrator)

    # Load state
    with open(os.path.join(restart_prefix, 'state.xml'), 'r') as f:
        state = XmlSerializer.deserialize(f.read())
    simulation.context.setState(state)

    # Selet atoms to save
    atom_indices = []
    top = md.load_pdb(pdbfile)
    res = [ r for r in top.topology.residues if r.name not in ("HOH", "NA", "CL") ]
    for r in res:
        for a in r.atoms:
            atom_indices.append(a.index)

    # Define reporter
    simulation.reporters.append(NetCDFReporter('traj.nc', netcdf_frequency, atomSubset=atom_indices))
    simulation.reporters.append(CheckpointReporter('checkpoint.chk', checkpoint_frequency))
    simulation.reporters.append(StateDataReporter('reporter.log', logging_frequency, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True))
    # Run
    simulation.step(nsteps)
    # Export state in xml format
    export_xml(simulation, system)
 

@click.command()
@click.option('--restart_prefix', default='../../../prep', help='path to load restart files')
@click.option('--initial_prefix', default='../../../prep', help='path to pdb file to load topology')
def cli(**kwargs):
    restart_prefix = kwargs["restart_prefix"]
    initial_prefix = kwargs["initial_prefix"]
    run(restart_prefix, initial_prefix)


if __name__ == "__main__":
    cli()