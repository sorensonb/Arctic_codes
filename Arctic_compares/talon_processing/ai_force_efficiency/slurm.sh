#!/usr/bin/env bash 
#-----------------------------------
# slurm-mm.sh
#
# David Apostal
# UND Computational Research Center
#
# Submit: 
#   sbatch slurm-mm.sh
# 
# Check status:
#   squeue | grep [userid]
#   squeue -u [userid]
#   squeue -j [jobid]
#-----------------------------------

#SBATCH --job-name=force_ai_test
#SBATCH --partition=gpu-code-test

# Sets the maximum time the job can run (hh:mm:ss).
#SBATCH --time=00:10:00
###SBATCH --time=01-12:00:00

# Specifies nodes for the job.
###SBATCH --nodes=1
### Slurm refers to cores as cpus. 
###SBATCH --cpus-per-task=16

# Sets the output filename.   %x = job name, %j = job-id
#SBATCH --output=%x.%j.txt

#SBATCH --mail-type=fail
#SBATCH --mail-user=blake.sorenson@und.edu

# This sets OMP_NUM_THREADS to the number of cores requested 
# from Slurm. OpenMP sets the number of threads per team 
# with OMP_NUM_THREADS. 
###export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "Number of OpenMP Threads: $OMP_NUM_THREADS"

##!#export OMP_PLACES=sockets
##!#export OMP_PROC_BIND=spread
##!#echo "OMP PLACES: $OMP_PLACES"
##!#echo "OMP PROC BIN: $OMP_PROC_BIND"

echo "Job start: $(date)"
echo ""

module load tensorflow2-py39-cuda11.8-gcc11/2.11.0


time python tensorflow_ai_test_noland103.py $1 $2 $3
#time python tensorflow_ai_test_noland102.py $1 $2 $3
#time python tensorflow_ai_test_noland101.py $1 $2 $3
#time python tensorflow_ai_test_noland100.py $1 $2 $3
#time python tensorflow_ai_test_noland99.py $1 $2 $3
#time python tensorflow_ai_test_noland98.py $1 $2 $3
#time python tensorflow_ai_test_noland97.py $1 $2 $3
#time python tensorflow_ai_test_noland96.py $1 $2 $3
#time python tensorflow_ai_test_noland95.py $1 $2 $3
#time python tensorflow_ai_test_noland94.py $1 $2 $3
#time python tensorflow_ai_test_noland93.py $1 $2 $3
#time python tensorflow_ai_test_noland92.py $1 $2 $3
#time python tensorflow_ai_test_noland91.py $1 $2 $3
#time python tensorflow_ai_test_noland90.py $1 $2 $3
#time python tensorflow_ai_test_noland89.py $1 $2 $3
#time python tensorflow_ai_test_noland88.py $1 $2 $3
#time python tensorflow_ai_test_noland86.py $1 $2 $3
#time python tensorflow_ai_test_noland85.py $1 $2 $3
#time python tensorflow_ai_test_noland84.py $1 $2 $3
#time python tensorflow_ai_test_noland83.py $1 $2 $3
#time python tensorflow_ai_test_noland82.py $1 $2 $3
#time python tensorflow_ai_test_noland81.py $1 $2 $3
#time python tensorflow_ai_test_noland80.py $1 $2 $3
#time python tensorflow_ai_test_noland79.py $1 $2 $3
#time python tensorflow_ai_test_noland78.py $1 $2 $3
#time python tensorflow_ai_test_noland77.py $1 $2 $3
#time python tensorflow_ai_test_noland76.py $1 $2 $3
#time python tensorflow_ai_test_noland75.py $1 $2 $3
#time python tensorflow_ai_test_noland74.py $1 $2 $3
#time python tensorflow_ai_test_noland73.py $1 $2 $3
#time python tensorflow_ai_test_noland72.py $1 $2 $3
#time python tensorflow_ai_test_noland71.py $1 $2 $3
#time python tensorflow_ai_test_noland70.py $1 $2 $3
#time python tensorflow_ai_test_noland69.py $1 $2 $3
#time python tensorflow_ai_test_noland68.py $1 $2 $3
#time python tensorflow_ai_test_noland67.py $1 $2 $3
#time python tensorflow_ai_test_noland66.py $1 $2 $3
#time python tensorflow_ai_test_noland65.py $1 $2 $3
#time python tensorflow_ai_test_noland64.py $1 $2 $3
#time python tensorflow_ai_test_noland63.py $1 $2 $3
#time python tensorflow_ai_test_noland62.py $1 $2 $3
#time python tensorflow_ai_test_noland61.py $1 $2 $3
#time python tensorflow_ai_test_noland60.py $1 $2 $3
#time python tensorflow_ai_test_noland59.py $1 $2 $3
#time python tensorflow_ai_test_noland58.py $1 $2 $3
#time python tensorflow_ai_test_noland57.py $1 $2 $3
#time python tensorflow_ai_test_noland56.py $1 $2 $3
#time python tensorflow_ai_test_noland55.py $1 $2 $3
#time python tensorflow_ai_test_noland54.py $1 $2 $3
#time python tensorflow_ai_test_noland53.py $1 $2 $3
#time python tensorflow_ai_test_noland52.py $1 $2 $3
#time python tensorflow_ai_test_noland51.py $1 $2 $3
#time python tensorflow_ai_test_noland50.py $1 $2 $3
#python model_test.py
#time python tensorflow_ai_test.py $1 $2
#time python tensorflow_ai_test_noland49.py $1 $2 $3

echo ""
echo "Job end: $(date)"

