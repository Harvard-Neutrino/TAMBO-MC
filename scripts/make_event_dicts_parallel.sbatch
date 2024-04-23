#!/bin/bash
#SBATCH -J TAMBO_trigger_sim
#SBATCH -c 1 # Number of cores
#SBATCH -p shared,arguelles_delgado
#SBATCH --mem 8000
#SBATCH -t 0-01:30:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/holyscratch01/arguelles_delgado_lab/Everyone/jlazar/TAMBO_event_dicts/%A_%a.out # Standard output
#SBATCH -e /n/holyscratch01/arguelles_delgado_lab/Everyone/jlazar/TAMBO_event_dicts/%A_%a.err # Standard output
#SBATCH --array 1-49,51-400

export PATH=$PATH:/n/home12/jlazar/julias/julia-1.9.2/bin/
export JULIA_DEPOT_PATH=/n/home12/jlazar/julias/julia-1.9.2

LENGTH=2000
DELTAS=150
printf -v PADDED_IDX "%05d" $SLURM_ARRAY_TASK_ID
BASEDIR="/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jlazar/graphnet_simulation/00000/00000_${PADDED_IDX}"
OUTFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jlazar/graphnet_simulation/00000/00000_${PADDED_IDX}/event_dicts_00000_${PADDED_IDX}.jld2"

CMD="time julia /n/home12/jlazar/TAMBO-MC/scripts/make_event_dicts_parallel.jl --length ${LENGTH} --deltas ${DELTAS} --outfile ${OUTFILE} --basedir ${BASEDIR}"
echo $CMD
$CMD
