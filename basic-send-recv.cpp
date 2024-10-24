#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <mpi.h>

typedef enum {
    TAG_WORKER_DISPATCH,
    TAG_WORKER_RETURN,
} MpiTags;

int do_work(int init_rank) {
    int i;
    for (i=0; i<1000; i++) ;
    return init_rank + i;
}

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize "work" variables
    int* ret_vals = (int*)malloc(size * sizeof(int));
    int work_val = -1;

    //Set up precision timer
    using std::chrono::milliseconds;
    auto t1 = std::chrono::high_resolution_clock::now();

    /**
     * The master process sends each worker an index (address) into the ret_vals array, initialized with worker rank:
     * worker_0:ret_vals[1] --> worker_1:work_val
     * worker_0:ret_vals[2] --> worker_2:work_val
     * ...
     * worker_0:ret_vals[N] --> worker_N:work_val
     */
    if (!rank) { // rank-0 is master
        for (int dest=1; dest < size; ++dest) { // 1, ... size-1
            ret_vals[dest] = dest; // initialize: "you are number X"
            MPI_Send(&(ret_vals[dest]), 1, MPI_INT, dest, TAG_WORKER_DISPATCH, MPI_COMM_WORLD);
        }
        // Master process "initializes" its work_val equal to its rank
        work_val = 0;
    } else {
        // "Initialize" work_val to whatever is received from master process
        MPI_Recv(&work_val, 1, MPI_INT, 0, TAG_WORKER_DISPATCH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    /**
     * Processes each "do some work"...
     * None of the worker processes are "aware" that work_val contains their rank
     */
    work_val = do_work(work_val);

    /**
     * Workers write their individual work_val back to master's ret_vals array at their respective indices:
     * worker_1:work_val --> worker_0:ret_vals[1]
     * worker_2:work_val --> worker_0:ret_vals[2]
     * ...
     * worker_N:work_val --> worker_0:ret_vals[N]
     */
    if (!rank) { // rank-0 is master
        for (int source=1; source < size; ++source) { // 1, ... size-1
            MPI_Recv(&(ret_vals[source]), 1, MPI_INT, source, TAG_WORKER_RETURN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // Master process performs a different operation than the workers, none of which are "aware" of this array
        ret_vals[0] = work_val;
    } else {
        // Send the result of "doing some work" back to master process
        MPI_Send(&work_val, 1, MPI_INT, 0, TAG_WORKER_RETURN, MPI_COMM_WORLD);
    }

    // is this required?
    MPI_Barrier(MPI_COMM_WORLD);

    // Stop the clock!
    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed_time_ms = std::chrono::duration<double, std::milli> (t2 - t1).count();

    // Print the "results"
    if (!rank) {
        for (int worker_rank=0; worker_rank < size; ++worker_rank) {
            printf("%d ", ret_vals[worker_rank]);
        }
        printf("\n");
    }

    // Clean up
    free(ret_vals);

    // Finalize MPI environment
    MPI_Finalize();

    // Print net runtime of running <STUFF>
    if(!rank) {
        std::cout << std::left<< std::setfill(' ')
                  << " size= " << std::setw(12) << std::setprecision(4) << size
                  << " elapsed= " << std::setw(12) << std::setprecision(4) << elapsed_time_ms << " ms "
                  << std::endl;
    }

    return 0;

}