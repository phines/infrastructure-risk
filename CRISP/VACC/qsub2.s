#!/bin/bash
#PBS -l nodes=1:eth10gb:ppn=8,pmem=2gb,pvmem=2gb
## running duration up to, shortq has a 3-hour limit
#PBS -l walltime=3:00:00
## queue to run in
#PBS -q workq
## Name of job
#PBS -N out_CSV
## Join STDERR TO STDOUT. (omit if you want separate STDOUT/STDERR output files)
#PBS -j oe
## Output to a specific file
#PBS -o /users/m/k/mkellygo/CRISP/output0/${PBS_JOBNAME}_${PBS_JOBID}.out
## Send mail on job start, end, aborted by batch system, or f for non-zero exit
#PBS -m bea
#PBS -M mkellygo@uvm.edu
#switcher julia = julia-1.1.0
ti=$(date +"%H:%M:%S.%N")
echo "started job at $ti"
cd /users/m/k/mkellygo/CRISP
module load julia
## times and runs program
time julia scripts/combine_resilience_csvs1.jl
ti=$(date +"%H:%M:%S.%N")
echo "completed job at ti"

