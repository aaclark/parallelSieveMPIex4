/* File:     seives.c
 * Purpose:  Find all the primes for the natural number in the range 1..Max
 *
 * Compile:  gcc seives.c -o seives -lpthread  -lm
 * Run:      ./seives <the highest ><max> <number of threads> <n>
 *           n is the number of terms of the series to use.
 *           n should be evenly divisible by the number of threads
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

int main (int argc, char ** argv) {
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
    gettimeofday(&start, NULL);
    bool* marked_natural_numbers = NULL;
    int* prime_numbers_sequential = NULL;
    int worker_from, worker_to, nr_of_prime_numbers_sequential, sequential_max;
    int* worker_ranges_from = (int*)malloc(size * sizeof(int));

    if (!rank) {
        // A boolean array where you can check if a natural number n (<=max) is 'marked' by inspecting marked_natural_numbers[n-1]
        // We want to mark all numbers that is not a prime
        marked_natural_numbers = (bool*)malloc(max * sizeof(bool));
        for (int i = 0; i < max; i++) {
            marked_natural_numbers[i] = false;
        }

        marked_natural_numbers[0] = true; // manually mark the first natural number 1
        int k = 2;
        sequential_max = ceil(sqrt(max)); // Round upward to ensure sequential_max*sequential_max <= max

        mark_all_numbers_not_prime_sequential(marked_natural_numbers, k, sequential_max);

        nr_of_prime_numbers_sequential = count_primes(marked_natural_numbers, sequential_max);

        printf("Number of primes found from the sequential part is %d \n", nr_of_prime_numbers_sequential);
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

            MPI_Isend(&worker_from, 1, MPI_INT, i, TAG_WORKER_FROM, MPI_COMM_WORLD, &request);
            MPI_Isend(&worker_to, 1, MPI_INT, i, TAG_WORKER_TO, MPI_COMM_WORLD, &request);
            MPI_Isend(&nr_of_prime_numbers_sequential, 1, MPI_INT, i, TAG_PRIME_COUNT, MPI_COMM_WORLD, &request);
            MPI_Isend(prime_numbers_sequential, nr_of_prime_numbers_sequential, MPI_INT, i, TAG_PRIMES, MPI_COMM_WORLD, &request);
            worker_from += chunk_size;
        }
    }

    //Let's update worker_from worker_to, prime_numbers_seqeuntial so each process can work on their range and mark the numbers that are not primes
    MPI_Recv(&worker_from, 1, MPI_INT, 0, TAG_WORKER_FROM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&worker_to, 1, MPI_INT, 0, TAG_WORKER_TO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&nr_of_prime_numbers_sequential, 1, MPI_INT, 0, TAG_PRIME_COUNT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    prime_numbers_sequential = (int*)malloc(nr_of_prime_numbers_sequential * sizeof(int));
    MPI_Recv(prime_numbers_sequential, nr_of_prime_numbers_sequential, MPI_INT, 0, TAG_PRIMES, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    printf("Process %d will work on chunk from %d to %d\n", rank, worker_from, worker_to);

    // Create a boolean subarray for the worker to mark the numbers that are not primes and send back the result to process 0
    // Process 0 can then match and update the sub array with the correct indexes at marked_natural_numbers
    int chunk_length = worker_to - worker_from + 1;
    bool* marked_natural_numbers_worker_chunk = (bool*)malloc(chunk_length * sizeof(bool));
    for (int j = 0; j < chunk_length; j++) {
        marked_natural_numbers_worker_chunk[j] = false; // Initialize to false
    }
    for (int natural_number = worker_from; natural_number <= worker_to; natural_number++) {
        for (int i = 0; i < nr_of_prime_numbers_sequential; i++) {
            int prime_number = prime_numbers_sequential[i];
            if (natural_number % prime_number == 0) {
                marked_natural_numbers_worker_chunk[natural_number-worker_from] = true;
                break;
            }
        }
    }

    printf("Process %d is starting to send it's computed result \n", rank);
    MPI_Request request;
    MPI_Isend(marked_natural_numbers_worker_chunk, chunk_length, MPI_C_BOOL, 0, TAG_CHUNK_RESULT, MPI_COMM_WORLD, &request);

    if (!rank) {
        int primes_pre = count_primes(marked_natural_numbers,max);
        for (int i = 0; i < size; i++) {
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, TAG_CHUNK_RESULT, MPI_COMM_WORLD, &status);
            int sender = status.MPI_SOURCE;
            int chunk_size_receive;
            int worker_range_from = worker_ranges_from[sender];
            if (sender==size-1) {
                chunk_size_receive=max-worker_ranges_from[sender]+1;
            } else {
                chunk_size_receive=worker_ranges_from[sender+1] - worker_range_from;
            }

            bool* marked_natural_numbers_chunk = (bool*)malloc(chunk_size_receive * sizeof(bool));
            MPI_Recv(marked_natural_numbers_chunk, chunk_size_receive, MPI_C_BOOL, sender, TAG_CHUNK_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(int i=0;i<chunk_size_receive;i++) {
                if (marked_natural_numbers_chunk[i]){
                    int marked_natural_number = i + worker_range_from;
                    marked_natural_numbers[marked_natural_number-1] = true;
                }
            }
            free(marked_natural_numbers_chunk);
            int primes_post = count_primes(marked_natural_numbers,max);
        }

        fflush(stdout);
        gettimeofday(&end, NULL);

        FILE* file_report = fopen("report.txt", "w");
        if (file_report == NULL) {
            printf("Error opening file for report!\n");
            return 1;
        }
        FILE* file_primes = fopen("primes.txt", "w");
        if (file_primes == NULL) {
            printf("Error opening file for writing the primes!\n");
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
        fprintf(file_report, "The program is using %d threads to find all the primes for the natural numbers in the range 1 to %d \n", size, max);
        fprintf(file_report, "Calculating all the primes took %f seconds:\n", total_time_spent);
        fprintf(file_report, "The sequential part calculated primes in the range 1 to %d and took %f seconds:\n", sequential_max, sequential_time_spent);
        fprintf(file_primes, "The first 1000 primes:\n");
        int number_of_primes = 0;
        for (int i = 0; i < max; i++) {
            // To avoid creating a huge file lets just write the first 1000 primes
            if (number_of_primes < 1000) {
                if (!marked_natural_numbers[i]){
                    fprintf(file_primes, "%d ", (i + 1));
                    number_of_primes++;
                }
            }
        }
        number_of_primes = count_primes(marked_natural_numbers, max);

        fclose(file_primes);
        fprintf(file_report, "Number of primes found %d\n", number_of_primes);
        fclose(file_report);
        fprintf(file_statistics, "\nUsing %d threads to find all the primes for the natural numbers in the range 1 to %d took %f seconds:\n", size, max, total_time_spent);
        fprintf(file_statistics, "The sequential part took %f seconds:\n", sequential_time_spent);
        fprintf(file_statistics, "The threaded part took %f seconds:\n", total_time_spent - sequential_time_spent);
        fclose(file_statistics);

        free(marked_natural_numbers);
        free(prime_numbers_sequential);
        printf("Number of primes found %d, more details can be read in report.txt \n", number_of_primes);
    }
    // Since we used non blocking send add a barrier to ensure that main thread has received
    // marked_natural_numbers_worker_chunk from all processes before we free it
    MPI_Barrier(MPI_COMM_WORLD);
    free(marked_natural_numbers_worker_chunk);

    MPI_Finalize();
    return 0;
}
