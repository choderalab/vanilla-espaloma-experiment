#!/bin/bash
#BSUB -P "benchmark"
#BSUB -J "@@@JOBNAME@@@"
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
#BSUB -W 12:00
###BSUB -m "ld-gpu ly-gpu lj-gpu ll-gpu lv-gpu"
#BSUB -o out_%J_%I.stdout
#BSUB -eo out_%J_%I.stderr
#BSUB -L /bin/bash

source ~/.bashrc
OPENMM_CPU_THREADS=1


# chnage dir
echo "changing directory to ${LS_SUBCWD}"
cd $LS_SUBCWD


# Report node in use
echo "======================"
hostname
env | sort | grep 'CUDA'
nvidia-smi
echo "======================"


# run job
conda activate perses-espaloma-0.3.0-v3
script_path=/home/takabak/data/espaloma-benchmark/benchmark-protein-ligand-vanilla/experiment/script
python ${script_path}/openmm_restart.py --restart_prefix @@@RESTART_PREFIX@@@
