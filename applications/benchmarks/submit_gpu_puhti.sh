#!/bin/bash -l
#SBATCH -J hila_bench
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-socket=2
#SBATCH --account=Project_2001973
#SBATCH --partition=gputest
#SBATCH --gres=gpu:v100:4
#

module load gcc/9.1.0 cuda/11.1.0 openmpi/4.0.5-cuda
#module load gcc/8.3.0 cuda/10.1.168 openmpi/4.0.3-cuda

UCX_MEMTYPE_CACHE=n srun bench_fermion > bench_fermion_4
UCX_MEMTYPE_CACHE=n srun bench_field > bench_field_4
UCX_MEMTYPE_CACHE=n srun bench_FFT > bench_FFT_4


