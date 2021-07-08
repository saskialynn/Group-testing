#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cdonnat@stanford.edu
#SBATCH --partition=hns,stat,normal,owners
#SBATCH --job-name=simulations
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --qos=normal
#SBATCH --time=48:00:00


ml load R/3.6.1 eigen gcc/6.3.0
cd $SCRATCH/Group_testing/Group-testing/

mode="uniform"
modeprev="none"
for tau in 0 #0.0001 0.001 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
       namefile="real_data_pool_$2_prev_$1_tau${tau}_${mode}_${modeprev}"
       Rscript computations.R $1 ${tau} $2 ${mode} ${modeprev} ${namefile}
done
