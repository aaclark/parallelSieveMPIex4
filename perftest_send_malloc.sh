#!/bin/zsh

# Output file to store results
output_file="benchmark_results_send_malloc.txt"

# Number of threads to test
proc_counts=(1 2 4 8 16 32)

# Max
max_n=(10000000 20000000 40000000 80000000 160000000 320000000)

# Clear the output file before starting the benchmarks
echo "" > $output_file

typeset -i N procs
# Run the benchmarks
for N in "${max_n[@]}"; do
  for procs in "${proc_counts[@]}"; do
    # Set the number of threads for OpenMP
    export OMP_NUM_THREADS=procs
    echo "P=${procs}, N=${N}"

    # Run the matrix_multiply program and capture the output
    result=$(mpiexec -np "$procs" ./sieve_send_malloc "$N")

    # Append the result to the output file
    echo "$result" >> $output_file
  done
  result=$(mpiexec -np 96 --hostfile hosts ./sieve_send_malloc "$N")
  echo "$result" >> $output_file
done

echo "Benchmarking complete! Results written to $output_file"
