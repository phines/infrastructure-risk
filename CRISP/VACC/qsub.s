#PBS -l nodes=1:ppn=6,pmem=4gb,pvmem=6gb
#PBS -l walltime=03:00:00
#PBS -q shortq
#PBS -joe
#PBS -N NEW_TEST_208
#switcher matlab = matlab2017b

time=$(date +"%H:%M:%S.%N")
echo "started job $PASS part $PASSB at $time"
cd $HOME/scratch/gauss_scrdir/dcsimsep-SHADE
# matlab -nodesktop -nodisplay -nosplash -nojvm -r "getRiskySupersetsVACC($PASS $PASSB); exit;"
matlab -nodesktop -nodisplay -nosplash -nojvm -r "SHADEdriver($PASSB); exit;"
time=$(date +"%H:%M:%S.%N")
echo "completed job $PASS part $PASSB at $time"
