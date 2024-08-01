#!/bin/bash
#
# PLD Slurm script
#
#SBATCH --job-name=PLD_jets
#SBATCH --partition=long
#SBATCH --mem=48G
#SBATCH --cpus-per-task=12
#SBATCH --time=10:00:00
#SBATCH --array=1-10
#
########################################################################
 
# Let's output which task we are within the job array
echo "First full run of PLD calculations on Slurm"
echo "Hello from Slurm job array task: $SLURM_ARRAY_TASK_ID"
 

echo "Checking memory"

# Get the job ID
JOB_ID=$SLURM_JOB_ID

# Query job details
MEM_LIMIT=$(scontrol show job "${JOB_ID}_${SLURM_ARRAY_TASK_ID}"| grep -oP 'mem=\K[^ ]+G')

echo "$MEM_LIMIT"
echo "$SLURM_ARRAY_TASK_ID"
echo "Beginning PLD calculation. the command being run is:"
echo "./PLDPermJob.py ${SLURM_ARRAY_TASK_ID} ${MEM_LIMIT}"
./PLDPermJob.py ${SLURM_ARRAY_TASK_ID} ${MEM_LIMIT}
echo "All done"