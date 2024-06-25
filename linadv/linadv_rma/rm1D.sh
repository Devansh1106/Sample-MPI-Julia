#!/bin/bash

sh -c "cat num_sol_par_*.txt > numerical_parallel.txt"
sh -c "cat exact_sol_par_*.txt > exact_parallel.txt"

sh -c "rm num_sol_par_*.txt"
sh -c "rm exact_sol_par_*.txt"