#!/bin/bash
#BSUB -P "analyze"
#BSUB -J "analyze"
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q cpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -W 0:10
#BSUB -o out_%J_%I.stdout
#BSUB -eo out_%J_%I.stderr
#BSUB -L /bin/bash

source ~/.bashrc
OPENMM_CPU_THREADS=1
#export OE_LICENSE=~/.openeye/oe_license.txt   # Open eye license activation/env

# change dir
echo "changing directory to ${LS_SUBCWD}"
cd $LS_SUBCWD

# Report node in use
echo "======================"
hostname
env | sort | grep 'CUDA'
nvidia-smi
echo "======================"

# conda
conda activate perses-espaloma-0.3.0-v3

# run
script_path=/home/takabak/data/espaloma-benchmark/benchmark-protein-ligand-vanilla/experiment/script
python ${script_path}/analysis.py
