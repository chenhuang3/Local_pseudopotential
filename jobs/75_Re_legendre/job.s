#!/bin/bash
#SBATCH -J "fhi"
#SBATCH -p huang_q
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH -o output.o
#SBATCH -e output.e

module purge 
module load intel openmpi

ulimit -s unlimited
cd $SLURM_SUBMIT_DIR/

./pswatch  k.ini  > log 

wait


