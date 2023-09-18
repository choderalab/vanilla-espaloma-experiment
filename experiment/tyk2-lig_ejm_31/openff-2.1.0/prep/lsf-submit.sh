#!/bin/bash
#BSUB -P "tyk2"
#BSUB -J "prep"
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
#BSUB -W  0:30
#BSUB -L /bin/bash
#BSUB -o out_%J_%I.stdout
#BSUB -eo out_%J_%I.stderr

source ~/.bashrc
OPENMM_CPU_THREADS=1

echo "changing directory to ${LS_SUBCWD}"
cd $LS_SUBCWD
conda activate perses-espaloma-0.3.0-v3

# Report node in use
hostname

# Report CUDA info
env | sort | grep 'CUDA'

# run
script_path=/home/takabak/data/espaloma-benchmark/benchmark-protein-ligand-vanilla/experiment/script
python ${script_path}/create_system.py --protein_file ../../crd/target.pdb --ligand_file ../../crd/ligands.sdf --small_molecule_forcefield openff-2.1.0.offxml
