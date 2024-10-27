/* File:     sieve_prime_count_reduce.c
 * Purpose:  Find all the primes for the natural number in the range 1..Max
 *
 * Compile:  mpicc sieves_prime_count.c -o sieves_prime_count -lm
 * Run:      mpiexec -np <n> ./seives <max>
 *           n is the number of processes
 *           We will find the primes in the range 1 to Max 
 * Output:   Stores the result in three files
 *           report.txt will report how many for primes was found and the execution times for this run
 *           primes.txt will contain all the primes from this run
 *           statistics text will append information for each run about execution times for the selected parameters
 */

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

#define MAX_NAME_SIZE 42

typedef enum {
  TAG_WORKER_FROM,
  TAG_WORKER_TO,
  TAG_PRIME_COUNT,
  TAG_PRIMES,
  TAG_CHUNK_RESULT
} MpiTags;

/* First a sequential solution for finding primes
   Mark all multiples of k between k*k and Max in the marked_natural_numbers boolean array.*/
void mark_all_multiples_of_k_sequential(bool marked_natural_numbers[], int k, int max) {
  int number_to_mark = k*k;
  while(number_to_mark <= max) {
    if (!marked_natural_numbers[number_to_mark - 1]) {
      marked_natural_numbers[number_to_mark-1]=true;
    }
    number_to_mark +=k;
  }
}

int find_next_smallest_unmarked_number_sequential(bool marked_natural_numbers[], int k, int max) {
  for (int i = k; i < max; i++) {
    if (!marked_natural_numbers[i - 1]) {
      return i;
    }
  }
    
  // It is assumed that when mark_all_numbers_not_prime_sequential calls this function a next unmarked number is always found.
  // So it shouldn't be possible to end up here but better to report and exit if we unexpectedly haven't exited the function yet.
  fprintf(stderr, "Error: Number not found.\n");
  exit(EXIT_FAILURE);
}
/*
  This function marks the natural numbers from k to max in the boolean array marked_natural_number that are not primes.
*/
void mark_all_numbers_not_prime_sequential(bool marked_natural_numbers[], int k, int max) {
  do {
    mark_all_multiples_of_k_sequential(marked_natural_numbers, k, max);
    k = find_next_smallest_unmarked_number_sequential(marked_natural_numbers, k+1, max);
  }  while (k*k <=max);
}

// Count how many primes there are (natural numbers that are not marked) in the range 1 .. max 
int count_primes(bool marked_natural_numbers[], int natural_number_max) {
  int primes_found = 0;
  for (int i = 0;i < natural_number_max; i++) {
    if (!marked_natural_numbers[i]) {
      primes_found++;
    }
  }
  return primes_found;
}

int find_nr_of_primes(int primes[], int primes_length, int from, int to) {
  int prime_count=0;
  for(int natural_number=from;natural_number<=to;natural_number++) { 
    bool is_prime=true;
    for (int i = 0; i < primes_length; i++) {
      int prime_number = primes[i];
      if (natural_number % prime_number == 0) {
	is_prime=false;
	break;
      }
    }
    if(is_prime) {
      prime_count++;
    }
  }
  return prime_count;
}

int main (int argc, char ** argv) {
  char name[MAX_NAME_SIZE];
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  int max = 110;
  if (argc > 1){
    max = atoi(argv[1]);
  }
  
  struct timeval start, starting_threads, end;
  bool* marked_natural_numbers = NULL;
  int* prime_numbers_sequential = NULL;
  int worker_from, worker_to, nr_of_prime_numbers_sequential;
  int* worker_ranges_from = (int*)malloc(size * sizeof(int));
  int total_prime_count=0;
  int sequential_max = ceil(sqrt(max)); // Round upward to ensure sequential_max*sequential_max <= max
  int total_range = max - sequential_max;
  int chunk_size = floor((total_range)/size);
  
  worker_from = sequential_max + 1;
  gettimeofday(&starting_threads, NULL);
  
  worker_from = sequential_max + 1 + chunk_size*rank;
  if (rank==size-1) {
    worker_to = max;
  } else {
    worker_to = worker_from + chunk_size - 1;
  }
  
  
  if (!rank) {
    gettimeofday(&start, NULL);
    // A boolean array where you can check if a natural number n (<=max) is 'marked' by inspecting marked_natural_numbers[n-1]
    
    // We want to mark all numbers that is not a prime
    marked_natural_numbers = (bool*)malloc(sequential_max * sizeof(bool));
    for (int i = 0; i < sequential_max; i++) {
      marked_natural_numbers[i] = false;
    }
    
    marked_natural_numbers[0] = true; // manually mark the first natural number 1
    int k = 2;
        
    mark_all_numbers_not_prime_sequential(marked_natural_numbers, k, sequential_max);

    nr_of_prime_numbers_sequential = count_primes(marked_natural_numbers, sequential_max);
    total_prime_count = nr_of_prime_numbers_sequential;
    printf("Number of primes found from the sequential part is %d\n", nr_of_prime_numbers_sequential);
    prime_numbers_sequential = (int*)malloc(nr_of_prime_numbers_sequential * sizeof(int));
    
    int prime_index = 0;
    for (int i=0; i<sequential_max; i++) {
      if (!marked_natural_numbers[i]) {
	int prime = i + 1;
	prime_numbers_sequential[prime_index] = prime;
	prime_index++;
      }
    }
  }

  gettimeofday(&starting_threads, NULL);
  // Step 1: Broadcast nr_of_prime_numbers_sequential
  MPI_Bcast(&nr_of_prime_numbers_sequential, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  // Step 2: Allocate prime_numbers_sequential on non-root ranks only
  if (rank != 0) {
    prime_numbers_sequential = (int*)malloc(nr_of_prime_numbers_sequential * sizeof(int));
  }
  
  // Step 3: Broadcast prime_numbers_sequential from rank 0 to all other ranks
  MPI_Bcast(prime_numbers_sequential, nr_of_prime_numbers_sequential, MPI_INT, 0, MPI_COMM_WORLD);
  
  printf("Process %d will work on chunk from %d to %d\n", rank, worker_from, worker_to);
  int chunk_prime_count = find_nr_of_primes(prime_numbers_sequential, nr_of_prime_numbers_sequential, worker_from, worker_to);
  MPI_Reduce(&chunk_prime_count, &total_prime_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   
  if (!rank) {
    total_prime_count+=nr_of_prime_numbers_sequential;
    
    fflush(stdout);
    gettimeofday(&end, NULL);
    
    FILE* file_report = fopen("report.txt", "w");
    if (file_report == NULL) {
      printf("Error opening file for report!\n");
      return 1;
    }
        
    FILE* file_statistics = fopen("statistics.txt", "a");
    if (file_statistics == NULL) {
      printf("Error opening file for statistics!\n");
      return 1;
    }

    fflush(stdout);
    double sequential_time_spent = ((double)(starting_threads.tv_usec - start.tv_usec) / 1000000.0) + starting_threads.tv_sec - start.tv_sec;
    double total_time_spent = ((double)(end.tv_usec - start.tv_usec) / 1000000.0) + end.tv_sec - start.tv_sec;
    fprintf(file_report, "The program is using %d processes to find all the primes for the natural numbers in the range 1 to %d \n", size, max);
    fprintf(file_report, "Calculating all the primes took %f seconds:\n", total_time_spent);
    fprintf(file_report, "The sequential part calculated primes in the range 1 to %d and took %f seconds:\n", sequential_max, sequential_time_spent);
    
    fprintf(file_report, "Number of primes found %d\n", total_prime_count);
    fclose(file_report);
    fprintf(file_statistics, "\nUsing %d processes to find all the primes for the natural numbers in the range 1 to %d took %f seconds:\n", size, max, total_time_spent);
    fprintf(file_statistics, "The sequential part took %f seconds:\n", sequential_time_spent);
    fprintf(file_statistics, "The parallel part took %f seconds:\n", total_time_spent - sequential_time_spent);
    fclose(file_statistics);
    
    free(marked_natural_numbers);
    free(prime_numbers_sequential);
    printf("Number of primes found %d, more details can be read in report.txt \n", total_prime_count);
  }
  // Since we used non blocking send add a barrier to ensure that main thread has received
  // marked_natural_numbers_worker_chunk from all processes before we free it 
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Finalize();
  return 0;
}
