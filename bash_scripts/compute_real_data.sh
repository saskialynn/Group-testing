#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cdonnat@stanford.edu
#SBATCH --partition=hns,stat,normal,owners
#SBATCH --job-name=simulations_graphs



cd $SCRATCH/Group_testing/experiments

for pool in {1..100}
do
for prev in 0.001 0.005 0.01 0.015 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.2 0.25 0.3
do
   #for tau in 0 0.0001 0.001 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
   #do
   #    namefile="real_data_pool_${pool}_prev_${prev}_tau${tau}"
   #    echo ${namefile}
       sbatch real_data.sh ${prev}  ${pool}  #${namefile}
  # done
done

done
