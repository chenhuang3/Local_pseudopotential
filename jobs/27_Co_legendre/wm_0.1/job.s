#!/bin/bash
#SBATCH --constraint='intel'
#SBATCH -J "psp"
#SBATCH -p backfill2
#SBATCH -n  1
#SBATCH --time=4:00:00
#SBATCH -o output-global.o
#SBATCH -e output-global.e


module purge 
module load intel openmpi

ulimit -s unlimited
cd $SLURM_SUBMIT_DIR/

./pswatch Co.ini > log 

wait


