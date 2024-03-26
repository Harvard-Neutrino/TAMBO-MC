#!/bin/bash
#SBATCH -c 1 # Number of cores
#SBATCH -p arguelles_delgado,shared,arguelles_delgado_gpu
#SBATCH --mem 10000
#SBATCH -t 0-04:00 # Maximum execution time (D-HH:MM)
#SBATCH -J TAMBO_March26th_Test
#SBATCH -o /n/holyscratch01/arguelles_delgado_lab/Everyone/pzhelnin/corsika_tambo/March26th_GraphNet/output/output_%A.out
#SBATCH -e /n/holyscratch01/arguelles_delgado_lab/Everyone/pzhelnin/corsika_tambo/March26th_GraphNet/error/error_%A.err

export PATH=$PATH:/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/julia-1.9.1/bin
export PYTHONPATH="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/TauRunner:$PYTHONPATH"
singularity_command=$1
eval $singularity_command 