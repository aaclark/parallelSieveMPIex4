#!/bin/sh
mpicxx -o sieve_bcast ex1/sieve_bcast.cpp ;
mpiexec -n "$1" ./sieve_bcast "$2" # --hostfile hosts
