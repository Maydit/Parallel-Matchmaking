/*
* Parallel Matchmaking Fall 2019
* David May, Thomas Clarke, Kousuke Tachida
*/

//Includes

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<mpi.h>

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime
#endif

//defines
#define NUM_TICKS 1000 //Number of ticks
#ifdef F2
#define MM_Dist matchmaking_dist_f2
#else
#define MM_Dist matchmaking_dist_f1
#endif
typedef struct Player {
  int mmr;
  unsigned int ping;
  int true_mmr;
} Player;

//Global vars
unsigned int total_players = 1000; //1,000

//Function declarations
float matchmaking_dist_f1();
void input_number(int in);

//Main
int main(int argc, char ** argv) {
  int mpi_rank;
  int mpi_size;
  double start_time = 0;
  double setup_time = 0;
  double total_time = 0;
  //Init MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  //Handle argvs & IO
  FILE * outfile = NULL;

  if(mpi_rank == 0) {
    outfile = fopen("data.txt", "w");
  }

  //Get time
  if(mpi_rank == 0) {
    start_time = GetTimeBase();
  }

  //Alloc memory

  //Report setup time
  if(mpi_rank == 0) {
    setup_time = GetTimeBase() - start_time;
    fprintf(outfile, "Total setup time:\t%012.6f\n", setup_time);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(int i = 0; i < NUM_TICKS; ++i) {
    //Exchange data

    //Run matchmaking algorithm

    //Report statistics
    MPI_Barrier(MPI_COMM_WORLD);
  }
  //Get total time
  if(mpi_rank == 0) {
    total_time = GetTimeBase() - start_time;
    fprintf(outfile, "Total runtime:\t\t%012.6f\nAverage time per tick:\t%012.6f\n", total_time, (total_time - setup_time) / NUM_TICKS);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  //close files
  if(mpi_rank == 0) {
    fclose(outfile);
  }
  //free memory
  return 0;
}

//Other functions here

//distance function
float matchmaking_dist_f1() {
  return 0.0;
}

void input_number(int in) {
  switch(in) {
    case 1:
      total_players = 10000; //10,000
      break;
    case 2:
      total_players = 100000; //100,000
      break;
    case 3:
      total_players = 1000000; //1,000,000
      break;
  }
}