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
#define GAME_LENGTH 5 //Number of ticks per game

#ifdef ELO
#define MM_Dist matchmaking_dist_elo
#define START_MMR 1200
#define nextmmr next_mmr_elo
#else
#define MM_Dist matchmaking_dist_trueskill
#define START_MMR 2500
#define nextmmr next_mmr_trueskill
#endif

#define START_UNCERT 830

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct Player {
  int mmr;
  unsigned int uncertainty;
  unsigned int ping;
  int true_mmr;
  float lenience;
  int playing;
  int group_size;
  int group_leader;
} Player;

//Global vars
unsigned int total_players = 1000; //1,000

//Function declarations
float matchmaking_dist_elo();
float matchmaking_dist_trueskill();
void input_number(int in);
int next_mmr_elo();
int next_mmr_trueskill();
unsigned int nextping();
float playerDistance(Player p1, Player p2);

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


  // setup
  float threshold = 1.0;
  total_players = 1024;
  int MAX_RANK_PLAYERS = total_players / mpi_size;
  int N_EXCHANGE_PLAYERS = MAX_RANK_PLAYERS/8;
  int N_RANK_PLAYERS = 6*N_EXCHANGE_PLAYERS;
  int GAME_PLAYERS = 2;
  int GAME_LENGTH = 5;
  int N_EXCHANGED;

  //Alloc memory
  Player * player_arr = malloc(sizeof(Player) * MAX_RANK_PLAYERS);
  //init
  srand48(mpi_rank);
  for(int i = 0; i < MAX_RANK_PLAYERS; ++i) {
    player_arr[i].ping        = nextping();
    player_arr[i].mmr         = START_MMR;
    player_arr[i].uncertainty = START_UNCERT;
    player_arr[i].true_mmr    = nextmmr();
    player_arr[i].playing     = 0; //not playing
    player_arr[i].lenience    = 1.0; //initial lenience
  }
  // initialize player distance
  float **distance = malloc(sizeof(float*)*MAX_RANK_PLAYERS);
  for (int i = 0; i < MAX_RANK_PLAYERS; i++) {
    distance[i] = malloc(sizeof(float)*MAX_RANK_PLAYERS);
  }
  for (int p = 0; p < MAX_RANK_PLAYERS; p++) {
    for (int j = p+1; j < MAX_RANK_PLAYERS; j++) {
      distance[p][j] = playerDistance(player_arr[p], player_arr[j]);
    }
  }
  
  //Report setup time
  if(mpi_rank == 0) {
    setup_time = GetTimeBase() - start_time;
    fprintf(outfile, "Total setup time:\t%012.6f\n", setup_time);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(int i = 0; i < NUM_TICKS; ++i) {
    // try to match players not in a match
    for (int p = 0; p < MAX_RANK_PLAYERS-1; p++) {
      if (player_arr[p].playing > 0) {
        continue;
      }
      for (int j = p+1; j < MAX_RANK_PLAYERS; j++) {
        if (player_arr[j].playing > 0) {
          continue;
        }
        // if (match()) {
        //   // match players
        //   player_arr[p].group;
        // 
        // 
        //   // start game
        //   // player_arr[p].playing = GAME_LENGTH;
        //   // player_arr[j].playing = GAME_LENGTH;
        // }
      }
    }
    
    // for (int i = 0; i < MAX_RANK_PLAYERS; i++) {
    // 
    // }
    // increase threshold for each unmatched player
    

    // exchange players at lower bound of each rank

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
  free(player_arr);
  return 0;
}

//Other functions here

float playerDistance(Player p1, Player p2) {
  if (p1.mmr >= p2.mmr) {
    return (float) p1.mmr - p2.mmr;
  } else {
    return (float) p2.mmr - p1.mmr;
  }
}

//distance function
float matchmaking_dist_elo() {
  return 0.0;
}

float matchmaking_dist_mmr() {
  return 0.0;
}

int n2_cached = 0;
double n2 = 0.0;
float gaussian(float stdev, float mean) {
  //Normal distribution centered around 1200
  //With standard deviation 2000 / 7
  //https://stackoverflow.com/questions/19944111/creating-a-gaussian-random-generator-with-a-mean-and-standard-deviation
  // static double n2 = 0.0;
  // static int n2_cached = 0;
  if(!n2_cached) {
    double x, y, r;
    do {
      x = 2.0*drand48() - 1;
      y = 2.0*drand48() - 1;
      r = x*x + y*y;
    } while(r == 0 || r > 1.0);
    double d = sqrt(-2.0 * log(r) / r);
    double n1 = x*d;
    n2 = y*d;
    double res = mean + n1 * stdev;
    n2_cached = 1;
    return res;
  } else {
    n2_cached = 0;
    return mean + n2 * stdev;
  }
}

int next_mmr_elo() {
  return (int) gaussian(2000.0 / 7.0, START_MMR);
}

int next_mmr_trueskill() {
  return (int) gaussian(START_UNCERT, START_MMR);
}

unsigned int nextping() {
  float g = gaussian(25, 50.0);
  return MAX(1, (int) g);
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