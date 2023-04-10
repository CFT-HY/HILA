#!/bin/bash -l
#SBATCH --job-name=job1         # Job name
#SBATCH --output=out%j          # Name of stdout output file
#SBATCH --error=err%j           # Name of stderr error file
#SBATCH --partition=standard-g  # Partition (queue) name
#SBATCH --nodes=64              # Total number of nodes
#SBATCH --ntasks-per-node=8     # 8 MPI ranks per node, 8 total (2x8)
#SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
#SBATCH --exclude=nid005360,nid005359
#SBATCH --time=04:10:00         # Run time (d-hh:mm:ss)
#SBATCH --account=Project_XXX   # Project for billing

# if fftw is also needed
module load CrayEnv PrgEnv-cray craype-accel-amd-gfx90a cray-mpich rocm fftw

CPU_BIND="map_cpu:48,56,16,24,1,8,32,40"

export MPICH_GPU_SUPPORT_ENABLED=1


srun --cpu-bind=${CPU_BIND} ./program

