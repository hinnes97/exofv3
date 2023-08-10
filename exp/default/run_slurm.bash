#!/bin/bash
#SBATCH --ntasks=24
#SBATCH -o test_default.out
#SBATCH -e test_defualt.err
#SBATCH -w atmnode018
#SBATCH --partition=priority-rp 
run.bash
