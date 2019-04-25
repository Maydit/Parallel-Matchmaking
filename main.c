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
#define GAME_LENGTH 10 //Number of ticks per game

#ifdef ELO
#define MM_update matchmaking_update_elo
#define START_MMR 1500
#define BUCKET_SIZE 400
#define N_BUCKETS 7
#define nextmmr next_mmr_elo
#else
#define MM_update matchmaking_update_trueskill
#define START_MMR 2500
#define BUCKET_SIZE 500
#define nextmmr next_mmr_trueskill
#endif

#define START_UNCERT 830
#define MM_CONST 1.0f / 10.0f
#define PING_CONST 1.0f / 50.0f

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MPI_TAG_ANY 1

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

//Global vars
unsigned int total_players = 4096; //1,024

//Function declarations
void input_number(int in);
int next_mmr_elo();
int next_mmr_trueskill();
unsigned int nextping();
float playerDistance(Player p1, Player p2);
int matchmaking_update_elo(Player * a, Player * b);
void matchmaking_update_trueskill(Player * a, Player * b);
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
  int MAX_RANK_PLAYERS = total_players / mpi_size;
  int PCT_EXCH = 4;
  int N_EXCHANGE_PLAYERS = MAX_RANK_PLAYERS / PCT_EXCH;

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
    player_arr[i].lenience    = 1; //initial lenience
    player_arr[i].wait_time   = -1; //not playing
    player_arr[i].total_wait_time = 0; // total timesteps waited
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

  // file setup
  // int a_won;
  // FILE *fp_games;
  // MPI_File mpi_f_games;
  // char file[100];
  // int DATA_SIZE = 17;
  // int game_data[DATA_LENGTH];
  // 
  // char expNum[3] = "0";
  // if (mpi_rank == 0) {
  //   strcpy(file, "out/games_");
  //   strcat(file, expNum);
  //   strcat(file, ".dat");
  //   fp_games = fopen(file, "w");
  // }
  // strcpy(file, "out/games_");
  // strcat(file, expNum);
  // strcat(file, ".dat");
  // MPI_Status status;
  // MPI_File_delete(file, MPI_INFO_NULL);
  // MPI_Barrier(MPI_COMM_WORLD);
  // MPI_File_open(MPI_COMM_WORLD, file, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_f_games);
  // int MAX_GAMES = (MAX_RANK_PLAYERS/2)*((float)NUM_TICKS/(float)GAME_LENGTH);
  // printf("%d\n", max_games);
  // int n_games = 0;
  // MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Offset rank_offset = sizeof(int)*(mpi_rank*MAX_GAMES*DATA_LENGTH);
  // for (int i = 0; i < 4; i++) {
  //   // game i
  //   game_data[0] = i;
  //   game_data[1] = i*2;
  // 
  // }
  

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
            // n_games++;
            player_arr[j].wait_time = GAME_LENGTH;
            player_arr[j].lenience  = 1;
            player_arr[k].wait_time = GAME_LENGTH;
            player_arr[k].lenience  = 1;
            // record_data(game_data, &player_arr[j], &player_arr[k], i);
            
            MM_update(&player_arr[j], &player_arr[k]);
            
            
            // MPI_File_write_at(mpi_f_games, rank_offset+(n_games*DATA_LENGTH*sizeof(int)), game_data, data_size, MPI_INT, &status);
            
            
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
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if(mpi_rank == 0) {
    sort_players(player_arr, 0, MAX_RANK_PLAYERS - 1);
    for(int i = 0; i < MAX_RANK_PLAYERS; ++i) {
      pretty_print_debug_player(outfile, player_arr[i], i);
    }
  }
  
  int USE_REDUCE = 0;
  // calculate statistics
  float all_average_wait_time[mpi_size];
  float average_wait_time;
  int total_wait_time = 0;
  int total_games = 0;
  
  int r_i; // all ranks index
  int b; // mmr_bucket
  
  int sum_counts[N_BUCKETS];
  int counts[N_BUCKETS];
  int all_counts[N_BUCKETS*mpi_size];
  
  float average_wait_times[N_BUCKETS]; // holds result
  float wait_times[N_BUCKETS]; // gather source
  float all_wait_times[N_BUCKETS*mpi_size]; // gather dest
  
  float average_winrates[N_BUCKETS];
  float winrates[N_BUCKETS]; // gather source
  float all_winrates[N_BUCKETS*mpi_size]; // gather dest
  
  float average_mmrs[N_BUCKETS];
  float mmrs[N_BUCKETS];
  float all_mmrs[N_BUCKETS*mpi_size];
  
  float average_true_mmrs[N_BUCKETS];
  float true_mmrs[N_BUCKETS];
  float all_true_mmrs[N_BUCKETS*mpi_size];
  
  // initialize values as 0
  for (int b = 0; b < N_BUCKETS; b++) {
    counts[b] = 0;
    wait_times[b] = 0;
    winrates[b] = 0;
    mmrs[b] = 0;
    true_mmrs[b] = 0;
  }
  
  for (int i = 0; i < MAX_RANK_PLAYERS; i++) {
    // total_wait_time += player_arr[i].total_wait_time;
    // total_games += player_arr[i].games;
  
    b = get_bucket(player_arr[i].mmr);
    counts[b] += 1;
    wait_times[b] += (float)player_arr[i].total_wait_time/(float)player_arr[i].games;
    winrates[b] += (float)player_arr[i].wins/(float)player_arr[i].games;
    mmrs[b] += (float)player_arr[i].mmr;
    true_mmrs[b] += (float)player_arr[i].true_mmr;
  }
  
  // average_wait_time = (float)total_wait_time/(float)total_games;
  // MPI_Gather(&average_wait_time, 1, MPI_FLOAT, all_average_wait_time, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  if (USE_REDUCE) {
    for (int b = 0; b < N_BUCKETS; b++) {
      MPI_Reduce(&counts[b], &sum_counts[b], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&wait_times[b], &average_wait_times[b], 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&winrates[b], &average_winrates[b], 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&mmrs[b], &average_mmrs[b], 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&true_mmrs[b], &average_true_mmrs[b], 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  } else {
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
  
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  float bucket_counts;
  if (mpi_rank == 0) {
    // for (int i = 0; i < mpi_size; i++) {
    //   printf("%d %f\n", mpi_rank, all_average_wait_time[i]);
    // }
    for (int b = 0; b < N_BUCKETS; b++) {
      if (USE_REDUCE) {
        bucket_counts = (float)sum_counts[b];
        average_wait_times[b] /= bucket_counts;
        average_winrates[b] /= bucket_counts;
        average_mmrs[b] /= bucket_counts;
        average_true_mmrs[b] /= bucket_counts;
        printf("%d %d\t%.2f\t%.2f\t%.2f\t%.2f\n", b, sum_counts[b], average_wait_times[b], average_winrates[b], average_mmrs[b], average_true_mmrs[b]);
      } else {
        // initialize values
        counts[b] = 0;
        wait_times[b] = 0;
        winrates[b] = 0;
        mmrs[b] = 0;
        true_mmrs[b] = 0;
        // sum over the ranks
        for (int r = 0; r < mpi_size; r++) {
          r_i = r*N_BUCKETS + b;
          if (all_counts[r_i] > 0) {
            counts[b] += all_counts[r_i];
            wait_times[b] += all_wait_times[r_i];
            winrates[b] += all_winrates[r_i];
            mmrs[b] += all_mmrs[r_i];
            true_mmrs[b] += all_true_mmrs[r_i];
          }
        }
        bucket_counts = (float)counts[b];
        wait_times[b] /= bucket_counts;
        winrates[b] /= bucket_counts;
        mmrs[b] /= bucket_counts;
        true_mmrs[b] /= bucket_counts;
        printf("%d %d\t%.2f\t%.2f\t%.2f\t%.2f\n", b, counts[b], wait_times[b], winrates[b], mmrs[b], true_mmrs[b]);
      }
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
  return (int) gaussian(2000.0 / 4.0, START_MMR);
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
      total_players *= 16; //16,384
      break;
    case 2:
      total_players *= 32; //32,768
      break;
    case 3:
      total_players *= 64; //65,536
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
    return ((a.lenience + b.lenience)/2.0 > (MM_CONST * abs(a.mmr - b.mmr) + PING_CONST * (a.ping + b.ping)));
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

void matchmaking_update_trueskill(Player * a, Player * b) {
  //TODO
}

void pretty_print_debug_player(FILE * f, Player a, int i) {
  fprintf(f, "Player %d stats: \n\tMMR: \t\t%d\n\tTRUE MMR:\t%d\n\tPING:\t\t%d\n\tWAIT TIME:\t%d\n\tWINS:\t\t%d\n\tGAMES:\t\t%d\n\tWR:\t\t%d\n", i, a.mmr, a.true_mmr, a.ping, a.wait_time, a.wins, a.games, (int) ((float) a.wins / (float) a.games * 100));
}


void write_end_data(Player* players) {
  
}