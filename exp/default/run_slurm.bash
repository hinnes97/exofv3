#!/bin/bash
#SBATCH --n 24
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH -w atmnode018
#SBATCH --partition=priority-rp 

source user_config

START=0
END=8
INC=2

for day in $(seq $START $INC $((END-INC)) ) ; do
    run.bash -s -d $INC -r $day -i $net2/output_newfms/$run_name/$day/RESTART
done

