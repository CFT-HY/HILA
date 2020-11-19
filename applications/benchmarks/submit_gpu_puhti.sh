#!/bin/bash -l
#SBATCH -J su2adj
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --account=Project_2001973
#SBATCH --partition=gputest
#SBATCH --gres=gpu:v100:2
#

module load gcc/9.1.0 cuda/11.1.0 openmpi/4.0.5-cuda
#module load gcc/8.3.0 cuda/10.1.168 openmpi/4.0.3-cuda

UCX_MEMTYPE_CACHE=n srun bench_fermion > bench_fermion_2
UCX_MEMTYPE_CACHE=n srun bench_field > bench_field_2


