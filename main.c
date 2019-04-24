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

#define ELO

//defines
#define NUM_TICKS 10000 //Number of ticks
#define GAME_LENGTH 20 //Number of ticks per game

#ifdef ELO
#define MM_update matchmaking_update_elo
#define START_MMR 1000
#define nextmmr next_mmr_elo
#else
#define MM_update matchmaking_update_trueskill
#define START_MMR 2500
#define nextmmr next_mmr_trueskill
#endif

#define START_UNCERT 830
#define MM_CONST 1.0f / 10.0f
#define PING_CONST 1.0f / 10.0f

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MPI_TAG_ANY 1

typedef struct Player {
  int mmr;
  int uncertainty;
  int ping;
  int true_mmr;
  int wait_time; //time remaining
  int wins; //games won
  int games; //games played
  int lenience;
} Player;

//Global vars
unsigned int total_players = 16; //1,024

//Function declarations
void input_number(int in);
int next_mmr_elo();
int next_mmr_trueskill();
unsigned int nextping();
float playerDistance(Player p1, Player p2);
void matchmaking_update_elo(Player * a, Player * b);
void matchmaking_update_trueskill(Player * a, Player * b);
void sort_players(Player * arr, int low, int high);
int matchable(Player a, Player b);

void pretty_print_debug_player(FILE * f, Player a, int i);

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
  total_players = 1024;
  int MAX_RANK_PLAYERS = total_players / mpi_size;
  int N_EXCHANGE_PLAYERS = MAX_RANK_PLAYERS / 8;

  MPI_Datatype mpi_player_type;
  MPI_Type_contiguous(8, MPI_INT, &mpi_player_type);
  MPI_Type_commit(&mpi_player_type);

  //Alloc memory
  Player * player_arr = malloc(sizeof(Player) * MAX_RANK_PLAYERS);
  //init
  srand48(mpi_rank);
  for(int i = 0; i < MAX_RANK_PLAYERS; ++i) {
    player_arr[i].ping        = nextping();
    player_arr[i].mmr         = START_MMR;
    player_arr[i].uncertainty = START_UNCERT;
    player_arr[i].true_mmr    = nextmmr();
    player_arr[i].lenience    = 1; //initial lenience
    player_arr[i].wait_time   = -1; //not playing
    player_arr[i].games       = 0;
    player_arr[i].wins        = 0;
  }
  Player * ghost_arr_prev = malloc(sizeof(Player) * N_EXCHANGE_PLAYERS);
  Player * ghost_arr_next = malloc(sizeof(Player) * N_EXCHANGE_PLAYERS);
  
  MPI_Request * reqs_next_r = malloc(sizeof(MPI_Request) * N_EXCHANGE_PLAYERS);
  MPI_Request * reqs_prev_r = malloc(sizeof(MPI_Request) * N_EXCHANGE_PLAYERS);
  MPI_Request * reqs_next_s = malloc(sizeof(MPI_Request) * N_EXCHANGE_PLAYERS);
  MPI_Request * reqs_prev_s = malloc(sizeof(MPI_Request) * N_EXCHANGE_PLAYERS);
  MPI_Status  * stat_next   = malloc(sizeof(MPI_Status)  * N_EXCHANGE_PLAYERS);
  MPI_Status  * stat_prev   = malloc(sizeof(MPI_Status)  * N_EXCHANGE_PLAYERS);

  //Report setup time
  if(mpi_rank == 0) {
    setup_time = GetTimeBase() - start_time;
    fprintf(outfile, "Total setup time:\t%012.6f\n", setup_time);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(int i = 0; i < NUM_TICKS; ++i) {
    //sort
    sort_players(player_arr, 0, MAX_RANK_PLAYERS - 1);
    //exchange
    //send N_EXCHANGE_PLAYERS to prev
    if(mpi_rank > 0) {
      MPI_Isend(player_arr, N_EXCHANGE_PLAYERS, mpi_player_type, mpi_rank - 1, MPI_TAG_ANY, MPI_COMM_WORLD, &reqs_next_s[mpi_rank]); //send to prev
      MPI_Irecv(ghost_arr_next, N_EXCHANGE_PLAYERS, mpi_player_type, mpi_rank - 1, MPI_TAG_ANY, MPI_COMM_WORLD, &reqs_prev_r[mpi_rank]); //recv from prev
      MPI_Wait(&reqs_prev_r[mpi_rank], &stat_prev[mpi_rank]);
      for(int j = 0; j < N_EXCHANGE_PLAYERS; ++j) {
        player_arr[7 * N_EXCHANGE_PLAYERS + j] = ghost_arr_next[j];
      }
    }
    //send N_EX_PL to next
    if(mpi_rank < mpi_size - 1) {
      MPI_Isend(&player_arr[7 * N_EXCHANGE_PLAYERS], N_EXCHANGE_PLAYERS, mpi_player_type, mpi_rank + 1, MPI_TAG_ANY, MPI_COMM_WORLD, &reqs_prev_s[mpi_rank]); //send to next
      MPI_Irecv(ghost_arr_prev, N_EXCHANGE_PLAYERS, mpi_player_type, mpi_rank + 1, MPI_TAG_ANY, MPI_COMM_WORLD, &reqs_next_r[mpi_rank]); //recieve from next
      MPI_Wait(&reqs_next_r[mpi_rank], &stat_next[mpi_rank]);
      for(int j = 0; j < N_EXCHANGE_PLAYERS; ++j) {
        player_arr[j] = ghost_arr_prev[j];
      }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    // try to match players not in a match
    for(int j = 0; j < MAX_RANK_PLAYERS; ++j) {
      if(player_arr[j].wait_time < 0) {
        for(int k = j + 1; k < MAX_RANK_PLAYERS; ++k) {
          if(matchable(player_arr[j], player_arr[k])) {
            player_arr[j].wait_time = GAME_LENGTH;
            player_arr[k].wait_time = GAME_LENGTH;
            MM_update(&player_arr[j], &player_arr[k]);
          }
        }
      }
    }
    for(int j = 0; j < MAX_RANK_PLAYERS; ++j) {
      if(player_arr[j].wait_time >= 0) {
        player_arr[j].wait_time -= 1;
      } else {
        player_arr[j].lenience += 1;
      }
    }
    //Report statistics
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if(mpi_rank == 0) {
    sort_players(player_arr, 0, MAX_RANK_PLAYERS - 1);
    for(int i = 0; i < MAX_RANK_PLAYERS; ++i) {
      pretty_print_debug_player(outfile, player_arr[i], i);
    }
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

//Quicksort
void swap_h(Player * a, Player * b) {
  Player t = *a;
  *a = *b;
  *b = t;
}

int partition_h(Player * arr, int low, int high) {
  int pivot = arr[high].mmr;
  int i = (low - 1);
  for(int j = low; j <= high - 1; j++) {
    if(arr[j].mmr <= pivot) {
      i++;
      swap_h(&arr[i], &arr[j]);
    }
  }
  swap_h(&arr[i+1], &arr[high]);
  return (i+1);
}

void sort_players(Player * arr, int low, int high) {
  if(low < high) {
    int pi = partition_h(arr, low, high);
    sort_players(arr, low, pi - 1);
    sort_players(arr, pi + 1, high);
  }
}

int matchable(Player a, Player b) {
  if(a.wait_time < 0 && b.wait_time < 0) {
    //printf("leniences: %d, %d\n", a.lenience, b.lenience);
    if(a.lenience + b.lenience > MM_CONST * abs(a.mmr - b.mmr) + PING_CONST * (a.ping + b.ping)) return 1;
  }
  return 0;
}

void matchmaking_update_elo(Player * a, Player * b) {
  //get result of match
  float qa = pow(10, a->true_mmr / 400.0f);
  float qb = pow(10, b->true_mmr / 400.0f);
  float ea = qa / (qa + qb);
  float eb = qb / (qa + qb);
  int res_a = drand48() < ea;
  int res_b = !res_a;
  //update
  qa = pow(10, a->mmr / 400.0f);
  qb = pow(10, b->mmr / 400.0f);
  ea = qa / (qa + qb);
  eb = qb / (qa + qb);
  a->mmr = a->mmr + 32.0 * (res_a - ea);
  b->mmr = b->mmr + 32.0 * (res_b - eb);
  //
  a->wins += res_a;
  b->wins += res_b;
  a->games += 1;
  b->games += 1;
}

void matchmaking_update_trueskill(Player * a, Player * b) {
  //TODO
}

void pretty_print_debug_player(FILE * f, Player a, int i) {
  fprintf(f, "Player %d stats: \n\tMMR: \t\t%d\n\tTRUE MMR:\t%d\n\tPING:\t\t%d\n\tWAIT TIME:\t%d\n\tWINS:\t\t%d\n\tGAMES:\t\t%d\n\tWR:\t\t%d\n", i, a.mmr, a.true_mmr, a.ping, a.wait_time, a.wins, a.games, (int) ((float) a.wins / (float) a.games * 100));
}