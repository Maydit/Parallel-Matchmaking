# Parallel-Matchmaking
RPI Parallel Programming Spring 2019 Final Project

Use make to build on local and make BGQ to build on blue gene

Edit run.sh as you please

(Temporary link to writeup)
https://www.overleaf.com/project/5cb8f8a6c16dba589acd63ad

TODOS:

Run code to produce graphs

Write results

Write abstract

Write conclusion

Finish writing algo if not enough space

Make the code pretty (and test again)

///BELOW HERE this is the actual readme we will be submitting

Contributers:
David May, Thomas Clarke, Kousuke Tachida

David May contributed general algorithm along with mmr systems implementation. He also wrote sections of the paper on the mmr systems and the algorithm.

Thomas Clarke contributed research into algorithms and contributed majorly to the research paper.

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
