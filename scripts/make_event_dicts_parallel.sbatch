#!/bin/bash
#SBATCH -J TAMBO_trigger_sim
#SBATCH -c 1 # Number of cores
#SBATCH -p shared
#SBATCH --mem 8000
#SBATCH -t 0-05:00:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/holyscratch01/arguelles_delgado_lab/Everyone/jlazar/TAMBO_event_dicts/%A_%a.out # Standard output
#SBATCH -e /n/holyscratch01/arguelles_delgado_lab/Everyone/jlazar/TAMBO_event_dicts/%A_%a.err # Standard output
#SBATCH --array 26-50

export PATH=$PATH:/n/home12/jlazar/julia-1.9.2/bin/
export JULIA_DEPOT_PATH=/n/home12/jlazar/julia-1.9.2

CMD="time julia /n/home12/jlazar/TAMBO-MC/scripts/make_event_dicts_parallel.jl --length ${LENGTH} --deltas ${DELTAS} --outfile ${OUTFILE} --basedir ${BASEDIR} --nparallel 50 --njob ${SLURM_ARRAY_TASK_ID} --run_desc ${RUN_DESC}"
echo $CMD
$CMD
