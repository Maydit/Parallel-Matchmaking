#!/bin/sh
module load xl
srun -N1 --ntasks-per-node=64 --overcommit -o -output.log ./main.xl
#call with $ sbatch --partition {PARTITION} --nodes {TOTAL NODES} --time {TIME} run.sh