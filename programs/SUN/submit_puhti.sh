#!/bin/bash
#SBATCH --job-name=SUN
#SBATCH --account=Project_2001973
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1

module load pgi

time srun SUN_GPU > output_GPU
time srun SUN_CPU > output_CPU

