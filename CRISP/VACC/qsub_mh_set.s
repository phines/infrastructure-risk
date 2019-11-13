#!/bin/bash
#PBS -l nodes=1:eth10gb:ppn=8,pmem=6gb,pvmem=6gb
## running duration up to, shortq has a 3-hour limit
#PBS -l walltime=3:00:00
## queue to run in
#PBS -q workq
## Name of job
#PBS -N two-opt_test01
## Join STDERR TO STDOUT. (omit if you want separate STDOUT/STDERR output files)
#PBS -j oe
## Output to a specific file
#PBS -o /users/m/k/mkellygo/CRISP/output1/9/${PBS_JOBNAME}_${PBS_JOBID}.out
## Send mail on job start, end, aborted by batch system, or f for non-zero exit
## #PBS -m bea
## #PBS -M mkellygo@uvm.edu
#PBS -t 1-100
#switcher julia = julia-1.1.0
ti=$(date +"%H:%M:%S.%N")
echo "started job $PASSx part $PASSB at $ti"
cd /users/m/k/mkellygo/CRISP
module load julia
## times and runs program
time julia scripts/run_CRISP_Rdist_mh_set.jl ${PBS_ARRAYID} $PASSy &
time julia scripts/run_CRISP_Rdist_mh_set1.jl ${PBS_ARRAYID} $PASSy &
time julia scripts/run_CRISP_Rdist_mh_set2.jl ${PBS_ARRAYID} $PASSy &
time julia scripts/run_CRISP_Rdist_mh_set3.jl ${PBS_ARRAYID} $PASSy &
time julia scripts/run_CRISP_Rdist_mh_set4.jl ${PBS_ARRAYID} $PASSy &
time julia scripts/run_CRISP_Rdist_mh_set5.jl ${PBS_ARRAYID} $PASSy &
time julia scripts/run_CRISP_Rdist_mh_set6.jl ${PBS_ARRAYID} $PASSy &
time julia scripts/run_CRISP_Rdist_mh_set7.jl ${PBS_ARRAYID} $PASSy &
wait
echo all processes complete
ti=$(date +"%H:%M:%S.%N")
echo "completed job at ti" ## $PASSx part $PASSy at $ti"