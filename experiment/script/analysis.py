#!/usr/bin/env python
"""
"""
import os, sys, math
import numpy as np
import glob
import click
import mdtraj
import logging
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
import seaborn as sns
# logging
logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore")

# plot settings
pd.options.display.max_rows = None
pd.options.display.max_columns = None
pd.options.display.precision = 1
pd.options.display.float_format = '{:.1f}'.format
params = {'legend.fontsize': 24, 
          'font.size': 36, 
          'axes.labelsize': 36,
          'axes.titlesize': 36,
          'xtick.labelsize': 36,
          'ytick.labelsize': 36,
          'savefig.dpi': 300, 
          #'figure.figsize': [64, 8],
          'xtick.major.size': 10,
          'xtick.minor.size': 7,
          'ytick.major.size': 10,
          'ytick.minor.size': 7}
plt.rcParams.update(params)

# settings
UNIT_NM_TO_ANGSTROMS = 10
UNIT_PS_TO_NS = 1/1000
LOGGING_FREQUENCY = 1000   # Default: 1000 ps



def load_traj(input_prefix, idx):
    """
    Load trajectory
    """

    # load initial (reference) structure
    init_pdb = "{}/complex-solvated.pdb".format(input_prefix)
    init_traj = mdtraj.load(init_pdb)
    atom_indices = init_traj.topology.select('not (water or resname HOH or resname NA or resname CL)')
    init_traj = init_traj.atom_slice(atom_indices)
    init_traj.save_pdb('traj.pdb')

    # load trajectories
    n = len(glob.glob(f"../sim{idx}/md*/traj.nc"))
    ncfiles = [ os.path.join("..", "sim" + str(idx), "md" + str(i), "traj.nc") for i in range(1, n+1) ]
    traj = mdtraj.load(ncfiles, top=init_traj.topology)
    #traj.save_netcdf('traj.nc')

    # check number of atoms
    assert init_traj.topology.n_atoms == traj.topology.n_atoms

    # Recenter and apply periodic boundary conditions to the molecules in each frame of the trajectory
    protein_indices = init_traj.topology.select('protein')
    mols = init_traj.atom_slice(protein_indices).topology.find_molecules()
    traj = traj.image_molecules(anchor_molecules=mols, inplace=True)
    #traj.save_netcdf('traj_image.nc')

    # Align trajectory to protein CA atoms (exclude first and last 5 residues)
    protein_ca_indices = init_traj.topology.select('protein and name == CA')
    traj.superpose(init_traj, 0, atom_indices=protein_ca_indices[5:-5])
    traj.save_netcdf('traj.nc')

    return init_traj, traj


def compute_custom_rmsd(traj, init_traj, atom_indices):
    """
    Compute RMSD

    Note: 
        mdtraj.rmsd automatically realigns the trajectory based on given atom indices
        https://github.com/mdtraj/mdtraj/issues/948
    """
    
    # https://github.com/srnas/barnaba/blob/master/barnaba/functions.py
    r = np.sqrt(3*np.mean((traj.xyz[:, atom_indices, :] - init_traj.xyz[0, atom_indices, :])**2, axis=(1,2)))
    
    return r


def calc_rmsd(init_traj, traj):
    """
    Compute RMSD using custom function

    We will compute RMSD using a custom function because mdtraj.rmsd seems to automatically 
    realign the trajectory based on given atom indices (https://github.com/mdtraj/mdtraj/issues/948)

    Parameters
    ----------
    init_traj : mdtraj.trajectory
        Solvated complex structure before minimization
    traj : mdtraj.trajectory
        Trajectory pre-aligned to C-alpha atoms
    """

    protein_ca_indices = init_traj.topology.select('protein and name == CA')
    #rmsd_pro = mdtraj.rmsd(traj, init_traj, atom_indices=protein_ca_indices[5:-5])
    rmsd_pro = compute_custom_rmsd(traj, init_traj, protein_ca_indices[5:-5])
    rmsd_pro = np.array(rmsd_pro) * UNIT_NM_TO_ANGSTROMS
    np.savetxt('rmsd_protein.out', rmsd_pro, fmt='%.3f')

    ligand_indices = init_traj.topology.select('resname == MOL and not symbol H')
    #rmsd_lig = mdtraj.rmsd(traj, init_traj, atom_indices=ligand_indices)
    rmsd_lig = compute_custom_rmsd(traj, init_traj, ligand_indices)
    rmsd_lig = np.array(rmsd_lig) * UNIT_NM_TO_ANGSTROMS
    np.savetxt('rmsd_ligand.out', rmsd_lig, fmt='%.3f')

    # realign trajectory wrt binding pocket (4 angstroms from ligand)
    ligand_indices = init_traj.topology.select('resname == MOL and not symbol H')
    pocket_atom_indices = mdtraj.compute_neighbors(init_traj, 0.4, ligand_indices)[0]
    pocket_residue_indices = [ atom.residue.index for atom in init_traj.topology.atoms if atom.index in pocket_atom_indices ]
    pocket_residue_indices = set(pocket_residue_indices)
    pocket_bb_indices = []
    for residue in init_traj.topology.residues:
        if residue.index in pocket_residue_indices:
            for atom in residue.atoms:
                if atom.is_backbone:
                    pocket_bb_indices.append(atom.index)
    traj.superpose(init_traj, 0, atom_indices=pocket_bb_indices)
    #rmsd_lig_bb = mdtraj.rmsd(traj, init_traj, atom_indices=ligand_indices)
    rmsd_lig_bb = compute_custom_rmsd(traj, init_traj, ligand_indices)
    rmsd_lig_bb = np.array(rmsd_lig_bb) * UNIT_NM_TO_ANGSTROMS
    np.savetxt('rmsd_ligand_bb.out', rmsd_lig_bb, fmt='%.3f')

    return rmsd_pro, rmsd_lig_bb
    

def plot(rmsd_dict, type):
    """
    Plot RMSD
    """
    rmsd_max = max(rmsd_dict['sim1'][f'{type}'].max(), rmsd_dict['sim2'][f'{type}'].max())

    fig, ax = plt.subplots(figsize=(12, 6))
    for key in rmsd_dict.keys():
        rmsd = rmsd_dict[key]
        rmsd = rmsd[f'{type}']
        x = np.arange(1, len(rmsd)+1) * LOGGING_FREQUENCY * UNIT_PS_TO_NS
        # x-axis
        ax.set_xlabel(r'Time (ns)')
        #ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.set_xlim([0, x[-1]]) 
        # y-axis
        ax.set_ylabel(r'RMSD (${\rm \AA}$)')
        #ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25))
        ax.set_ylim([0, math.ceil(rmsd_max)])
        ax.plot(x, rmsd, lw=1, label=f'{key}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"rmsd_{type}.png")
    plt.savefig(f"rmsd_{type}.svg", dpi=2400)
    

def run(input_prefix):
    """
    Compute rmsd and plot
    """

    # check number of simulation trials
    n_trials = len(glob.glob("../sim*"))
    # compute rmsd for each simulation trial
    rmsd_dict = {}
    for idx in range(1, n_trials+1):
        init_traj, traj = load_traj(input_prefix, idx)
        rmsd_pro, rmsd_lig = calc_rmsd(init_traj, traj)
        rmsd_dict[f"sim{idx}"] = {"protein": rmsd_pro, "ligand": rmsd_lig}
    # plot
    plot(rmsd_dict, type="protein")
    plot(rmsd_dict, type="ligand")


@click.command()
@click.option('--input_prefix', default='../../prep', help='path to pdb file to load topology')
def cli(**kwargs):
    input_prefix = kwargs['input_prefix']
    run(input_prefix)


if __name__ == '__main__':
    cli()

