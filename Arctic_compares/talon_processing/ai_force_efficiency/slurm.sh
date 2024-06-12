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
#SBATCH --time=00:30:00

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


#python model_test.py
time python tensorflow_ai_test.py $1 $2

echo ""
echo "Job end: $(date)"

