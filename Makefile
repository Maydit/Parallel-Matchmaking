all: main.c
	mpicc -I. -Wall -O3 main.c -o main.out -lm

BGQ-ELO: main.c
	mpixlc -I. -Wall -O3 main.c -o main.xl -lm -D BGQ

BGQ-TS: main.c
	mpixlc -I. -Wall -O3 main.c -o main.xl -lm -D BGQ -D MMR