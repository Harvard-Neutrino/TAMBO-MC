#!/bin/bash
#SBATCH -J TAMBO_trigger_sim
#SBATCH -c 1 # Number of cores
#SBATCH -p shared,arguelles_delgado,sapphire 
#SBATCH --mem 18000
#SBATCH -t 0-15:30:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/pzhelnin/corsika_tambo/Larger_valley_event_dicts/small/output/%A_%a.out # Standard output
#SBATCH -e /n/netscratch/arguelles_delgado_lab/Everyone/pzhelnin/corsika_tambo/Larger_valley_event_dicts/small/error/%A_%a.err # Standard output
#SBATCH --array 1-50

export PATH=/n/home03/pzhelnin/.juliaup/bin/julia:$PATH

LENGTH=5000
DELTAS=150
ALTMIN=2300
ALTMAX=3700
SIZE="medium"

printf -v PADDED_IDX "%05d" $SLURM_ARRAY_TASK_ID
BASEDIR="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/hard_spectrum/shower/larger_valley_00000_${PADDED_IDX}"
OUTFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/hard_spectrum/event_dicts/retest_${SIZE}_larger_valley_event_dicts_00000_${PADDED_IDX}.jld2"
SIMFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/hard_spectrum/sim_files/larger_valley_00000_${PADDED_IDX}.jld2"
#BASEDIR="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/pavel_cross_check/Larger_Valley/hshower/larger_valley_00000_${PADDED_IDX}"
#SIMFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/pavel_cross_check/larger_valley_00000_${PADDED_IDX}.jld2"
#OUTFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/event_dicts/Will_${SIZE}_larger_valley_event_dicts_00000_${PADDED_IDX}.jld2"

CMD="time julia /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/TAMBO-MC/scripts/1_make_event_dicts/make_event_dicts_parallel.jl --length ${LENGTH} --deltas ${DELTAS} --outfile ${OUTFILE} --basedir ${BASEDIR} --altmin ${ALTMIN} --altmax ${ALTMAX} --simfile ${SIMFILE} --size ${SIZE}"
echo $CMD
echo $SLURM_ARRAY_TASK_ID
$CMD
