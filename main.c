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
#define TOTAL_PLAYERS 32768 //8192 //16384//1024 //1048576 //total players
#define NUM_TICKS 10000 //Number of ticks
#define GAME_LENGTH 10 //Number of ticks per game
#define USE_REDUCE 1
#define USE_GATHER 0

#define MMR
#ifndef MMR
#define MM_update matchmaking_update_elo
#define MM_CONST 1.0 / 50
#define START_MMR 1500
#define BUCKET_SIZE 400
#define N_BUCKETS 7
#define nextmmr next_mmr_elo
#else
#define MM_update matchmaking_update_trueskill
#define MM_CONST 1.0 / 100
#define START_MMR 2500
#define BUCKET_SIZE 500
#define N_BUCKETS 9
#define nextmmr next_mmr_trueskill
#endif

#define START_UNCERT 830
#define PING_CONST 1.0f / 200
#define BETA 4.0

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MPI_TAG_ANY 1

double g_processor_frequency = 1600000000.0; // processing speed for BG/Q

typedef struct Player {
  int mmr;
  int uncertainty;
  int ping;
  int true_mmr;
  int wait_time; //time remaining
  int total_wait_time;
  int wins; //games won
  int games; //games played
  int lenience;
} Player;

//Function declarations
void input_number(int in);
int next_mmr_elo();
int next_mmr_trueskill();
unsigned int nextping();
float playerDistance(Player p1, Player p2);
int matchmaking_update_elo(Player * a, Player * b);
int matchmaking_update_trueskill(Player * a, Player * b);
void sort_players(Player * arr, int low, int high);
int matchable(Player a, Player b);
void record_data(int* game_data, Player* a, Player* b);
void pretty_print_debug_player(FILE * f, Player a, int i);
int get_bucket(int mmr);

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
  int MAX_RANK_PLAYERS = TOTAL_PLAYERS / mpi_size;
  int PCT_EXCH = 4;
  int N_EXCHANGE_PLAYERS = MAX_RANK_PLAYERS / PCT_EXCH;
  int TICK_FREQ = NUM_TICKS/32;

  MPI_Datatype mpi_player_type;
  MPI_Type_contiguous(9, MPI_INT, &mpi_player_type);
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
    player_arr[i].lenience    = 5; //initial lenience
    player_arr[i].wait_time   = -1; //not playing
    player_arr[i].total_wait_time = 0; // total timesteps waited
    player_arr[i].games       = 0;
    player_arr[i].wins        = 0;
  }
  Player * ghost_arr_prev = malloc(sizeof(Player) * N_EXCHANGE_PLAYERS);
  Player * ghost_arr_next = malloc(sizeof(Player) * N_EXCHANGE_PLAYERS);
  
  MPI_Request * reqs_next_r = malloc(sizeof(MPI_Request) * mpi_size);
  MPI_Request * reqs_prev_r = malloc(sizeof(MPI_Request) * mpi_size);
  MPI_Request * reqs_next_s = malloc(sizeof(MPI_Request) * mpi_size);
  MPI_Request * reqs_prev_s = malloc(sizeof(MPI_Request) * mpi_size);
  MPI_Status  * stat_next   = malloc(sizeof(MPI_Status)  * mpi_size);
  MPI_Status  * stat_prev   = malloc(sizeof(MPI_Status)  * mpi_size);
  
  // setup statistics
  double stat_time;
  double stat_start_time;  
  float bucket_counts;
  int r_i; // all ranks index
  int b; // mmr_bucket
  
  int *sum_counts = malloc(sizeof(int)*N_BUCKETS);
  int *counts = malloc(sizeof(int)*N_BUCKETS);
  int *all_counts = malloc(sizeof(int)*N_BUCKETS*mpi_size);
  
  float *average_wait_times = malloc(sizeof(float)*N_BUCKETS);
  float *wait_times = malloc(sizeof(float)*N_BUCKETS);
  float *all_wait_times = malloc(sizeof(float)*N_BUCKETS*mpi_size);
  
  float *average_winrates = malloc(sizeof(float)*N_BUCKETS);
  float *winrates = malloc(sizeof(float)*N_BUCKETS);
  float *all_winrates = malloc(sizeof(float)*N_BUCKETS*mpi_size);
  
  float *average_mmrs = malloc(sizeof(float)*N_BUCKETS);
  float *mmrs= malloc(sizeof(float)*N_BUCKETS);
  float *all_mmrs = malloc(sizeof(float)*N_BUCKETS*mpi_size);
  
  float *average_true_mmrs = malloc(sizeof(float)*N_BUCKETS);
  float *true_mmrs = malloc(sizeof(float)*N_BUCKETS);
  float *all_true_mmrs = malloc(sizeof(float)*N_BUCKETS*mpi_size);
  
  //Report setup time
  if(mpi_rank == 0) {
    setup_time = (GetTimeBase() - start_time) / g_processor_frequency;
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
      MPI_Irecv(ghost_arr_prev, N_EXCHANGE_PLAYERS, mpi_player_type, mpi_rank - 1, MPI_TAG_ANY, MPI_COMM_WORLD, &reqs_prev_r[mpi_rank]); //recv from prev
      MPI_Wait(&reqs_prev_r[mpi_rank], &stat_prev[mpi_rank]);
    }
    //send N_EX_PL to next
    if(mpi_rank < mpi_size - 1) {
      MPI_Isend(&player_arr[(PCT_EXCH - 1) * N_EXCHANGE_PLAYERS], N_EXCHANGE_PLAYERS, mpi_player_type, mpi_rank + 1, MPI_TAG_ANY, MPI_COMM_WORLD, &reqs_prev_s[mpi_rank]); //send to next
      MPI_Irecv(ghost_arr_next, N_EXCHANGE_PLAYERS, mpi_player_type, mpi_rank + 1, MPI_TAG_ANY, MPI_COMM_WORLD, &reqs_next_r[mpi_rank]); //recieve from next
      MPI_Wait(&reqs_next_r[mpi_rank], &stat_next[mpi_rank]);
    }
    //Overwrite data
    if(mpi_rank > 0) {
      for(int j = 0; j < N_EXCHANGE_PLAYERS; ++j) {
        player_arr[j] = ghost_arr_prev[j];
      }
    }
    if(mpi_rank < mpi_size - 1) {
      for(int j = 0; j < N_EXCHANGE_PLAYERS; ++j) {
        player_arr[(PCT_EXCH - 1) * N_EXCHANGE_PLAYERS + j] = ghost_arr_next[j];
      }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    // try to match players not in a match
    for(int j = 0; j < MAX_RANK_PLAYERS; ++j) {
      if(player_arr[j].wait_time < 0) {
        for(int k = j + 1; k < MAX_RANK_PLAYERS; ++k) {
          if(matchable(player_arr[j], player_arr[k])) {
            player_arr[j].wait_time = GAME_LENGTH;
            player_arr[j].lenience  = 5;
            player_arr[k].wait_time = GAME_LENGTH;
            player_arr[k].lenience  = 5;
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
        player_arr[j].total_wait_time += 1;
      }
    }
    //Report statistics
    if (i % TICK_FREQ == 0) {
      if (mpi_rank == 0) {
        printf("\ntick %d\n", i);
      }
      // initialize values as 0
      for (int b = 0; b < N_BUCKETS; b++) {
        counts[b] = 0;
        wait_times[b] = 0;
        winrates[b] = 0;
        mmrs[b] = 0;
        true_mmrs[b] = 0;
      }
      // get stats for each bucket
      for (int i = 0; i < MAX_RANK_PLAYERS; i++) {
        b = get_bucket(player_arr[i].mmr);
        counts[b] += 1;
        wait_times[b] += (float)player_arr[i].total_wait_time/(float)player_arr[i].games;
        winrates[b] += (float)player_arr[i].wins/(float)player_arr[i].games;
        mmrs[b] += (float)player_arr[i].mmr;
        true_mmrs[b] += (float)player_arr[i].true_mmr;
      }
      if (USE_REDUCE) {
        MPI_Barrier(MPI_COMM_WORLD);
        if(mpi_rank == 0) {
          stat_start_time = GetTimeBase();
        }
        MPI_Reduce(counts, sum_counts, N_BUCKETS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(wait_times, average_wait_times, N_BUCKETS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(winrates, average_winrates, N_BUCKETS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(mmrs, average_mmrs, N_BUCKETS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(true_mmrs, average_true_mmrs, N_BUCKETS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (mpi_rank == 0) {
          for (int b = 0; b < N_BUCKETS; b++) {
            bucket_counts = (float)sum_counts[b];
            average_wait_times[b] /= bucket_counts;
            average_winrates[b] /= bucket_counts;
            average_mmrs[b] /= bucket_counts;
            average_true_mmrs[b] /= bucket_counts;
            printf("%d %d\t%.2f\t%.2f\t%.2f\t%.2f\n", b, sum_counts[b], average_wait_times[b], average_winrates[b], average_mmrs[b], average_true_mmrs[b]);
          }
          stat_time = (GetTimeBase() - stat_start_time) / g_processor_frequency;
          printf("reduce record time:\t%012.6f\n", stat_time);
        }
      }
      
      if (USE_GATHER) {
        MPI_Barrier(MPI_COMM_WORLD);
        if(mpi_rank == 0) {
          stat_start_time = GetTimeBase();
        }
        // gather counts
        MPI_Gather(counts, N_BUCKETS, MPI_INT, 
                   all_counts, N_BUCKETS, MPI_INT, 
                   0, MPI_COMM_WORLD);
        // gather 
        MPI_Gather(wait_times, N_BUCKETS, MPI_INT, 
                   all_wait_times, N_BUCKETS, MPI_INT, 
                   0, MPI_COMM_WORLD);
        // gather winrates into all_winrates
        MPI_Gather(winrates, N_BUCKETS, MPI_FLOAT, 
                   all_winrates, N_BUCKETS, MPI_FLOAT, 
                   0, MPI_COMM_WORLD);
        // gather mmrs into all_mmrs
        MPI_Gather(mmrs, N_BUCKETS, MPI_FLOAT, 
                   all_mmrs, N_BUCKETS, MPI_FLOAT, 
                   0, MPI_COMM_WORLD);               
        // gather true_mmrs into all_true_mmrs
        MPI_Gather(true_mmrs, N_BUCKETS, MPI_FLOAT, 
                   all_true_mmrs, N_BUCKETS, MPI_FLOAT, 
                   0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (mpi_rank == 0) {
          for (int b = 0; b < N_BUCKETS; b++) {
            // initialize values
            sum_counts[b] = 0;
            average_wait_times[b] = 0;
            average_winrates[b] = 0;
            average_mmrs[b] = 0;
            average_true_mmrs[b] = 0;
            // sum over the ranks
            for (int r = 0; r < mpi_size; r++) {
              r_i = r*N_BUCKETS + b;
              if (all_counts[r_i] > 0) {
                sum_counts[b] += all_counts[r_i];
                average_wait_times[b] += all_wait_times[r_i];
                average_winrates[b] += all_winrates[r_i];
                average_mmrs[b] += all_mmrs[r_i];
                average_true_mmrs[b] += all_true_mmrs[r_i];
              }
            }
            bucket_counts = (float)sum_counts[b];
            average_wait_times[b] /= bucket_counts;
            average_winrates[b] /= bucket_counts;
            average_mmrs[b] /= bucket_counts;
            average_true_mmrs[b] /= bucket_counts;
            printf("%d %d\t%.2f\t%.2f\t%.2f\t%.2f\n", b, sum_counts[b], average_wait_times[b], average_winrates[b], average_mmrs[b], average_true_mmrs[b]);
          }
          stat_time = GetTimeBase() - stat_start_time;
          printf("gather record time:\t%012.6f\n", stat_time);
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // record final stats
  if (mpi_rank == 0) {
    printf("\ntick %d\n", NUM_TICKS);
  }
  // initialize values as 0
  for (int b = 0; b < N_BUCKETS; b++) {
    counts[b] = 0;
    wait_times[b] = 0;
    winrates[b] = 0;
    mmrs[b] = 0;
    true_mmrs[b] = 0;
  }
  // get stats for each bucket
  for (int i = 0; i < MAX_RANK_PLAYERS; i++) {
    b = get_bucket(player_arr[i].mmr);
    counts[b] += 1;
    wait_times[b] += (float)player_arr[i].total_wait_time/(float)player_arr[i].games;
    winrates[b] += (float)player_arr[i].wins/(float)player_arr[i].games;
    mmrs[b] += (float)player_arr[i].mmr;
    true_mmrs[b] += (float)player_arr[i].true_mmr;
  }
  if (USE_REDUCE) {
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_rank == 0) {
      stat_start_time = GetTimeBase();
    }
    MPI_Reduce(counts, sum_counts, N_BUCKETS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(wait_times, average_wait_times, N_BUCKETS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(winrates, average_winrates, N_BUCKETS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(mmrs, average_mmrs, N_BUCKETS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(true_mmrs, average_true_mmrs, N_BUCKETS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      for (int b = 0; b < N_BUCKETS; b++) {
        bucket_counts = (float)sum_counts[b];
        average_wait_times[b] /= bucket_counts;
        average_winrates[b] /= bucket_counts;
        average_mmrs[b] /= bucket_counts;
        average_true_mmrs[b] /= bucket_counts;
        printf("%d %d\t%.2f\t%.2f\t%.2f\t%.2f\n", b, sum_counts[b], average_wait_times[b], average_winrates[b], average_mmrs[b], average_true_mmrs[b]);
      }
      stat_time = (GetTimeBase() - stat_start_time) / g_processor_frequency;
      printf("reduce record time:\t%012.6f\n", stat_time);
    }
  }
  
  if (USE_GATHER) {
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_rank == 0) {
      stat_start_time = GetTimeBase();
    }
    // gather counts
    MPI_Gather(counts, N_BUCKETS, MPI_INT, 
               all_counts, N_BUCKETS, MPI_INT, 
               0, MPI_COMM_WORLD);
    // gather 
    MPI_Gather(wait_times, N_BUCKETS, MPI_INT, 
               all_wait_times, N_BUCKETS, MPI_INT, 
               0, MPI_COMM_WORLD);
    // gather winrates into all_winrates
    MPI_Gather(winrates, N_BUCKETS, MPI_FLOAT, 
               all_winrates, N_BUCKETS, MPI_FLOAT, 
               0, MPI_COMM_WORLD);
    // gather mmrs into all_mmrs
    MPI_Gather(mmrs, N_BUCKETS, MPI_FLOAT, 
               all_mmrs, N_BUCKETS, MPI_FLOAT, 
               0, MPI_COMM_WORLD);               
    // gather true_mmrs into all_true_mmrs
    MPI_Gather(true_mmrs, N_BUCKETS, MPI_FLOAT, 
               all_true_mmrs, N_BUCKETS, MPI_FLOAT, 
               0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      for (int b = 0; b < N_BUCKETS; b++) {
        // initialize values
        sum_counts[b] = 0;
        average_wait_times[b] = 0;
        average_winrates[b] = 0;
        average_mmrs[b] = 0;
        average_true_mmrs[b] = 0;
        // sum over the ranks
        for (int r = 0; r < mpi_size; r++) {
          r_i = r*N_BUCKETS + b;
          if (all_counts[r_i] > 0) {
            sum_counts[b] += all_counts[r_i];
            average_wait_times[b] += all_wait_times[r_i];
            average_winrates[b] += all_winrates[r_i];
            average_mmrs[b] += all_mmrs[r_i];
            average_true_mmrs[b] += all_true_mmrs[r_i];
          }
        }
        bucket_counts = (float)sum_counts[b];
        average_wait_times[b] /= bucket_counts;
        average_winrates[b] /= bucket_counts;
        average_mmrs[b] /= bucket_counts;
        average_true_mmrs[b] /= bucket_counts;
        printf("%d %d\t%.2f\t%.2f\t%.2f\t%.2f\n", b, sum_counts[b], average_wait_times[b], average_winrates[b], average_mmrs[b], average_true_mmrs[b]);
      }
      stat_time = (GetTimeBase() - stat_start_time) / g_processor_frequency;
      printf("gather record time:\t%012.6f\n", stat_time);
    }
  }
  
  //Get total time
  if(mpi_rank == 0) {
    total_time = (GetTimeBase() - start_time) / g_processor_frequency;
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
  free(sum_counts);
  free(counts);
  free(all_counts);
  free(average_wait_times);
  free(wait_times);
  free(all_wait_times);
  free(average_winrates);
  free(winrates);
  free(all_winrates);
  free(average_mmrs);
  free(mmrs);
  free(all_mmrs);
  free(average_true_mmrs);
  free(true_mmrs);
  free(all_true_mmrs);
  return 0;
}

//Other functions here
int get_bucket(int mmr) {
  if (mmr < BUCKET_SIZE) {
    return 0;
  } else if (mmr >= BUCKET_SIZE*(N_BUCKETS-1)) {
    return N_BUCKETS-1;
  } else {
    return mmr / BUCKET_SIZE;
  }
}

int n2_cached = 0;
double n2 = 0.0;
float gaussian(float stdev, float mean) {
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
  return (int) gaussian(2000.0 / 4.0, START_MMR);
}

int next_mmr_trueskill() {
  return (int) gaussian(START_UNCERT, START_MMR);
}

unsigned int nextping() {
  float g = gaussian(25, 50.0);
  return MAX(1, (int) g);
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
    return ((a.lenience + b.lenience)/10.0 > (MM_CONST * abs(a.mmr - b.mmr) + PING_CONST * (a.ping + b.ping)));
  }
  return 0;
}

int matchmaking_update_elo(Player * a, Player * b) {
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
  return res_a;
}

float Psi(float x);
float Phi(float x);

float sigmoid(float x) {
  //1 = .8
  //-1 = .2
  //https://www.microsoft.com/en-us/research/project/trueskill-ranking-system/
  return 1.0 / (exp(-x * 1.3) + 1.0);
}

int matchmaking_update_trueskill(Player * a, Player * b) {
  float sig_a = a->uncertainty / 100.0f;
  float sig_b = b->uncertainty / 100.0f;
  float a_mmr = a->mmr / 100.0f;
  float b_mmr = b->mmr / 100.0f;
  float a_tr  = a->true_mmr / 100.0f;
  float b_tr  = b->true_mmr / 100.0f;
  float c = sqrt(2 * BETA * BETA + sig_a * sig_a + sig_b * sig_b);
  //calc winner
  //float draw = exp(-1 * (a_mmr - b_mmr) *(a_mmr - b_mmr) / 2 / c / c) * (sqrt(2 * BETA * BETA / c / c));
  float ea = sigmoid((a_tr - b_tr) / BETA);
  //printf("%f\n", ea);
  int res = drand48() < ea;
  //update
  if(res) {
    a->mmr = a->mmr + 100 * (sig_a * sig_a / c) * (Phi((a_mmr - b_mmr) / c));
    b->mmr = b->mmr - 100 * (sig_b * sig_b / c) * (Phi((a_mmr - b_mmr) / c));
    a->uncertainty = (int) (100 * sqrt(sig_a * (1 - sig_a * sig_a / c / c * Psi((a_mmr - b_mmr) / c))));
    b->uncertainty = (int) (100 * sqrt(sig_b * (1 - sig_b * sig_b / c / c * Psi((a_mmr - b_mmr) / c))));
  } else {
    a->mmr = a->mmr - 100 * (sig_a * sig_a / c) * (Phi((b_mmr - a_mmr) / c));
    b->mmr = b->mmr + 100 * (sig_b * sig_b / c) * (Phi((b_mmr - a_mmr) / c));
    a->uncertainty = (int) (100 * sqrt(sig_a * (1 - sig_a * sig_a / c / c * Psi((b_mmr - a_mmr) / c))));
    b->uncertainty = (int) (100 * sqrt(sig_b * (1 - sig_b * sig_b / c / c * Psi((b_mmr - a_mmr) / c))));
  }
  //
  a->wins += res;
  b->wins += !res;
  a->games += 1;
  b->games += 1;
  return res;
}

void pretty_print_debug_player(FILE * f, Player a, int i) {
  fprintf(f, "Player %d stats: \n\tMMR: \t\t%d\n\tTRUE MMR:\t%d\n\tPING:\t\t%d\n\tWAIT TIME:\t%d\n\tWINS:\t\t%d\n\tGAMES:\t\t%d\n\tWR:\t\t%d\n", i, a.mmr, a.true_mmr, a.ping, a.wait_time, a.wins, a.games, (int) ((float) a.wins / (float) a.games * 100));
}

float Phi(float x) {
  //Reverse engineered
  //https://www.moserware.com/assets/computing-your-skill/The%20Math%20Behind%20TrueSkill.pdf
  return MAX(1.0 - x, 0.1);
}

float Psi(float x) {
  //Reverse engineered
  //https://www.moserware.com/assets/computing-your-skill/The%20Math%20Behind%20TrueSkill.pdf
  return 1.0 - 1.0 / (exp(-2.0 * x + 2.0) + 1.0);
}