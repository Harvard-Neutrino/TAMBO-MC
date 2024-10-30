#!/bin/bash
#SBATCH -J TAMBO_make_event_dicts
#SBATCH -c 1 # Number of cores
#SBATCH -p arguelles_delgado,shared,sapphire
##SBATCH -p shared,sapphire
##SBATCH -p arguelles_delgado
#SBATCH --mem 8000
#SBATCH -t 0-04:00:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/holyscratch01/arguelles_delgado_lab/Everyone/thomwg11/TAMBO_event_dicts/%A_%a.out # Standard output
#SBATCH -e /n/holyscratch01/arguelles_delgado_lab/Everyone/thomwg11/TAMBO_event_dicts/%A_%a.err # Standard output
#SBATCH --array 1-100

#export PATH=$PATH:/n/home02/thomwg11/.juliaup/bin
#export JULIA_DEPOT_PATH=/n/home02/thomwg11/.julia

LENGTH=2000
DELTAS=100
ALTMIN=1850
ALTMAX=3850
printf -v PADDED_IDX "%05d" $SLURM_ARRAY_TASK_ID
#BASEDIR="/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jlazar/graphnet_simulation/00000/00000_${PADDED_IDX}"
#BASEDIR="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Jan7th2024_WhitePaper_300k_no_thin/"
BASEDIR="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/coordinate_tests/showers/"
#BASEDIR="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Oct16th2023_WhitePaper_300k/"
#BASEDIR="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Jan7th2024_WhitePaper_300k_no_thin/Jan23rd2024_WhitePaper_300k_no_thin_rerun/"
#OUTFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jlazar/graphnet_simulation/00000/00000_${PADDED_IDX}/event_dicts_00000_${PADDED_IDX}.jld2"
#OUTFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/thomwg11/graphnet_simulation_tests_jeff/event_dicts.jld2"
#OUTFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/Jan7th2024_WhitePaper_300k_no_thin/event_dicts/event_dicts.jld2"
OUTFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/coordinate_tests/event_dicts/event_dicts.jld2"
SIMFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/coordinate_tests/coordinate_tests_00000_00001.jld2"
#SIMFILE="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Oct16th2023_WhitePaper_300k.jld2"

#CMD="time julia /n/home02/thomwg11/tambo/TAMBO-MC/scripts/1_make_event_dicts/make_event_dicts_parallel.jl --length ${LENGTH} --deltas ${DELTAS} --altmin ${ALTMIN} --altmax ${ALTMAX} --outfile ${OUTFILE} --basedir ${BASEDIR} --simfile ${SIMFILE} --nparallel 100 --njob ${SLURM_ARRAY_TASK_ID}"
CMD="time julia /n/home02/thomwg11/tambo/TAMBO-MC/scripts/1_make_event_dicts/make_event_dicts_parallel.jl --config perfect_valley_2000_150.toml --nparallel 100 --njob ${SLURM_ARRAY_TASK_ID}"
echo $CMD
$CMD
