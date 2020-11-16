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

module load gcc/9.1.0 openmpi/4.0.5-cuda cuda/11.1.0

srun bench_fermion > bench_fermion_new

