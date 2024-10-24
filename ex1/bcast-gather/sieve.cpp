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

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <mpi.h>

#define MAX_N 1000000  // Define a maximum N value
#define WRITEOUT false // debugging...

typedef enum {
    TAG_RANGE,
    TAG_PRIMES,
    TAG_CHUNK_RESULT
} MpiTags;


// Iterate all primes within range
void mark_multiples_sequential(bool* marked_numbers, int max) {
    /**
     * Assume each integer p = 1...sqrt(MAX) is prime,
     * Then for each integer p = 1...sqrt(MAX):
     *      For each multiple i of p, i = p*p, p*p+p, p*p+2p, ... sqrt(MAX):
     *          Mark i as non-prime
     */
    for (int i = 2; i <= max; i++) marked_numbers[i] = true;  // Initialize all to true (prime)

    for (int p = 2; p <= max; p++) { // For each integer p...
        if (marked_numbers[p]) {    // IF it is prime...
            for (int i = p * p; i <= max; i += p) marked_numbers[i] = false;  // Mark multiples as non-prime
        }
    }
}


// Iterate over multiples of a given prime
void mark_multiples_parallel(bool *marked_numbers, int start, int prime, int limit) {
    /**
     * For each multiple i of p, i = p, 2p, 3p, ... limit:
     *      Mark i as non-prime
     */
    for (int i = start; i < limit; i += prime) {
        marked_numbers[i] = false;
    }
}


// Print help or such
void usage(char* argv[]) {
    std::cerr << "Usage: \n"
              << "mpiexec -n <processes>" << argv[0] << " <N> <=" << MAX_N << "\n"
              << argv[0] << " -h\n";
}


int main(int argc, char *argv[]) {
    int rank, size, N;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if ((argc != 2)||!strcmp(argv[1],"-h")) {
        usage(argv);
        MPI_Finalize();
        std::exit(EXIT_FAILURE);
    }
    N = atol(argv[1]);
    if (N > MAX_N) {
        if (rank == 0) std::cerr << "N is too large. Max allowed is "<< MAX_N << std::endl;
        MPI_Finalize();
        return 1;
    }
    int sqrt_N = (int)sqrt(N);

    // Each process has a unique rank and can compute its own start and end (inclusive) indices
    int range_start =   (rank)      * (N / size) + 2; // Offset
    int range_end =     (rank + 1)  * (N / size) + 1; // (rank + 1) * (N / size) + 2 - 1;
    if (rank == size - 1) range_end = N; // Truncate the last section

    // TODO do we need to allocate dynamically???
    int range_size = range_end - range_start + 1;
    bool prime_array[range_size]; // = malloc...

    // This array is UNIQUE TO EACH WORKER; initialize...
    for (int i = 0; i < range_size; i++) prime_array[i] = true;

    // Stack-allocate space for gathering results (only for rank 0)
    bool global_primes[MAX_N + 1];  // Large enough to hold results up to MAX_N
    if (!rank) { for (int i = 0; i <= N; i++) { global_primes[i] = false; }}  // Initialize global array


    //Set up precision timer
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    //Start the clock!
    auto t1 = high_resolution_clock::now();
    {
        // Stack-allocate array for small primes (handled by rank 0)
        bool small_primes[sqrt_N + 1];
        if (!rank) {
            mark_multiples_sequential(small_primes, sqrt_N);
        }

        // Broadcast small primes to all processes
        MPI_Bcast(small_primes, sqrt_N + 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        {
            /** PARALLEL SECTION */

            // Each process marks multiples of small primes in its range
            for (int p = 2; p <= sqrt_N; p++) { // Skip 0, 1
                if (small_primes[p]) {
                    // determine the smallest multiple of p that is >= range_start; and smallest_multiple at that offset for p
                    int smallest_multiple = (range_start % p == 0) ? range_start : range_start + (p - (range_start % p));
                    if (smallest_multiple == p) smallest_multiple += p; // Avoid marking the prime itself
                    mark_multiples_parallel(prime_array, smallest_multiple, p, range_end);
                }
            }

        }
        // Gather results back to the root process
        MPI_Gather(prime_array, range_size, MPI_CXX_BOOL, global_primes + range_start - 2, range_size, MPI_CXX_BOOL, 0,
                   MPI_COMM_WORLD);

    }
    // Stop the clock!
    auto t2 = high_resolution_clock::now();
    auto elapsed_time_ms = duration<double, std::milli> (t2 - t1).count();

    // Output the results (in rank 0)
    if (!rank) {
        for (int i = 2; i <= N; i++) {
            if (global_primes[i] && WRITEOUT) printf("%d ", i);
        }
        printf("\n");
    }

    // Finalize the MPI environment
    MPI_Finalize();

    if(!rank) {
        // Print results and net runtime of running <STUFF>
        std::cout << std::left<< std::setfill(' ')
                  << " N= " << std::setw(12) << std::setprecision(4) << N
                  << " elapsed= " << std::setw(12) << std::setprecision(4) << elapsed_time_ms << " ms "
                  << std::endl;
    }

    return 0;
}
