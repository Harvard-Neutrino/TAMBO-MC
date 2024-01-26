#!/bin/bash
#SBATCH -J TAMBO_trigger_sim
#SBATCH -c 1 # Number of cores
#SBATCH -p shared
#SBATCH --mem 8000
#SBATCH -t 0-05:00:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/holyscratch01/arguelles_delgado_lab/Everyone/jlazar/TAMBO_trigger_sim/%A_%a.out # Standard output
#SBATCH -e /n/holyscratch01/arguelles_delgado_lab/Everyone/jlazar/TAMBO_trigger_sim/%A_%a.err # Standard output
#SBATCH --array 1-17

export PATH=$PATH:/n/home12/jlazar/julia-1.9.2/bin/
export JULIA_DEPOT_PATH=/n/home12/jlazar/julia-1.9.2
CMD="time julia /n/home12/jlazar/TAMBO-MC/scripts/trigger_sim_parallel.jl --length ${LENGTH} --deltas ${DELTAS} --outfile ${OUTFILE} --basedir ${BASEDIR} --nparallel 17 --njob ${SLURM_ARRAY_TASK_ID}"
echo $CMD
$CMD
