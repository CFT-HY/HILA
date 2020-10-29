#!/bin/bash -l
#SBATCH -J su2adj
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 00:10:00
#SBATCH -n 2
#SBATCH --account=Project_2001973
#SBATCH --partition=gpu 
#SBATCH --gres=gpu:v100:2
#

module load gcc/8.3.0 hpcx-mpi/2.5.0-cuda cuda/10.1.168

mpirun -n 2 bench_fermion

