#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cdonnat@stanford.edu
#SBATCH --partition=hns,stat,normal,owners
#SBATCH --job-name=simulations
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --qos=normal
#SBATCH --time=24:00:00


ml load R/3.6.1 eigen gcc/6.3.0
cd $SCRATCH/Group_testing/experiments

Rscript compute_p.R $1 $2 $3
