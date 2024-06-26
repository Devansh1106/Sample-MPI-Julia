#!/bin/bash

# Concatenate solution files horizontly 
sh -c "paste -d '\t\t' num_sol2D_par_*.txt > _numerical_parallel.txt"
sh -c "paste -d '\t\t' exact_sol2D_par_*.txt > _exact_parallel.txt" 

# Remove solution files generated by each rank
sh -c "rm num_sol2D_par_*.txt"
sh -c "rm exact_sol2D_par_*.txt"

# Concatenate mesh files in y direction vertically
sh -c "cat mesh_y_*.txt > mesh_y.txt"
sh -c "rm mesh_y_*.txt"

# Concatenate all files horizontly to a single file
sh -c "paste -d '\t\t' mesh_x.txt mesh_y.txt _numerical_parallel.txt > numerical_parallel.txt"
sh -c "paste -d '\t\t' mesh_x.txt mesh_y.txt _exact_parallel.txt > exact_parallel.txt"

# Removing unecessary files
sh -c "rm _numerical_parallel.txt"
sh -c "rm _exact_parallel.txt"
sh -c "rm mesh_*.txt"