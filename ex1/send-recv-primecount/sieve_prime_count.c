/* File:     seives.c
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
    // printf("Testing %d\n", natural_number);
    bool is_prime=true;
    for (int i = 0; i < primes_length; i++) {
      int prime_number = primes[i];
      if (natural_number % prime_number == 0) {
	is_prime=false;
	break;
      }
    }
    if(is_prime) {
      // printf("Prime found %d\n", natural_number);
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
  
  struct timeval start, starting_threads, end, gather_start, gather_end;
  bool* marked_natural_numbers = NULL;
  int* prime_numbers_sequential = NULL;
  int worker_from, worker_to, nr_of_prime_numbers_sequential, sequential_max;
  int* worker_ranges_from = (int*)malloc(size * sizeof(int));
  int local_prime_count=0;//Each process counts its primes
  int total_prime_count;
  if (!rank) {
    gettimeofday(&start, NULL);
    // A boolean array where you can check if a natural number n (<=max) is 'marked' by inspecting marked_natural_numbers[n-1]
    sequential_max = ceil(sqrt(max)); // Round upward to ensure sequential_max*sequential_max <= max
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
    
    int total_range = max - sequential_max;
    int chunk_size = floor((total_range)/size);
    worker_from = sequential_max + 1;
    gettimeofday(&starting_threads, NULL);
    MPI_Request request;
    for (int i = 0; i < size; i++) {
      worker_ranges_from[i] = worker_from;
      if (i==size-1) {
        worker_to = max;
      } else {
	worker_to = worker_from + chunk_size - 1;
      }

      // Send the message to itself non blocking to avoid deadlock
      if (i==rank) {
	MPI_Isend(&worker_from, 1, MPI_INT, i, TAG_WORKER_FROM, MPI_COMM_WORLD, &request);
	MPI_Isend(&worker_to, 1, MPI_INT, i, TAG_WORKER_TO, MPI_COMM_WORLD, &request);
	MPI_Isend(&nr_of_prime_numbers_sequential, 1, MPI_INT, i, TAG_PRIME_COUNT, MPI_COMM_WORLD, &request);
	MPI_Isend(prime_numbers_sequential, nr_of_prime_numbers_sequential, MPI_INT, i, TAG_PRIMES, MPI_COMM_WORLD, &request);
      } else {
	MPI_Send(&worker_from, 1, MPI_INT, i, TAG_WORKER_FROM, MPI_COMM_WORLD);
	MPI_Send(&worker_to, 1, MPI_INT, i, TAG_WORKER_TO, MPI_COMM_WORLD);
	MPI_Send(&nr_of_prime_numbers_sequential, 1, MPI_INT, i, TAG_PRIME_COUNT, MPI_COMM_WORLD);
	MPI_Send(prime_numbers_sequential, nr_of_prime_numbers_sequential, MPI_INT, i, TAG_PRIMES, MPI_COMM_WORLD);
      }
      worker_from += chunk_size;
    }
  }

  // Let's update worker_from worker_to, prime_numbers_seqeuntial so each process can work on their range and mark the numbers that are not primes
  MPI_Recv(&worker_from, 1, MPI_INT, 0, TAG_WORKER_FROM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&worker_to, 1, MPI_INT, 0, TAG_WORKER_TO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&nr_of_prime_numbers_sequential, 1, MPI_INT, 0, TAG_PRIME_COUNT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  prime_numbers_sequential = (int*)malloc(nr_of_prime_numbers_sequential * sizeof(int));
  MPI_Recv(prime_numbers_sequential, nr_of_prime_numbers_sequential, MPI_INT, 0, TAG_PRIMES, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  printf("Process %d will work on chunk from %d to %d\n", rank, worker_from, worker_to);
  MPI_Barrier(MPI_COMM_WORLD);
  int chunk_prime_count = find_nr_of_primes(prime_numbers_sequential, nr_of_prime_numbers_sequential, worker_from, worker_to);
  int len;
  MPI_Get_processor_name(name, &len);
  printf("Process %d running at host %s is starting to send it's computed result \n", rank, name);
  MPI_Request request;
  MPI_Isend(&chunk_prime_count, 1, MPI_INT, 0, TAG_CHUNK_RESULT, MPI_COMM_WORLD, &request);
  
  if (!rank) {
    for (int i = 0; i < size; i++) {
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, TAG_CHUNK_RESULT, MPI_COMM_WORLD, &status);      
      int sender = status.MPI_SOURCE;
      int nr_of_primes_found_by_process;
      MPI_Recv(&nr_of_primes_found_by_process, 1, MPI_INT, sender, TAG_CHUNK_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      total_prime_count += nr_of_primes_found_by_process;
    }
    
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