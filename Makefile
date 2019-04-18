all: main.c
	mpicc -I. -Wall -O3 main.c -o main.out

BGQ: main.c
	mpicc -I. -Wall -O5 main.c -o main.xl -D BGQ