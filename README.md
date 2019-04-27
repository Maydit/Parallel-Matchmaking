# Parallel-Matchmaking
RPI Parallel Programming Spring 2019 Final Project

[Formal writeup here](./writeup.pdf)

Contributers:
David May, Thomas Clarke, Kousuke Tachida

David May contributed general algorithm along with mmr systems implementation. He also wrote sections of the paper on the mmr systems and the algorithm.

Thomas Clarke contributed research into algorithms, contributed majorly to the research paper, and ran majority of BGQ tests.

Kousuke Tachida contributed to the statistics collection along with creating the graphs for the paper.

How to run:

scp main.c and Makefile to Blue Gene / Q

$ module load xl

$ make BGQ-ELO

OR

$ make BGQ-TS

depending on which one you want to test

$ emacs run.sh

write the following into run.sh:

'''

#!/bin/sh

module load xl

srun -N1 --ntasks-per-node=64 --overcommit -o -output.log ./main.xl

'''

Call the following command:

$sbatch --partition small --nodes 1 --time 25 run.sh

Here, if you want to change the size, change nodes 1 and N1 to some multiple of 2.

Output is stored in the slurm.log and the data.txt
