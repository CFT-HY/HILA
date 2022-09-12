#!/bin/bash

#SBATCH --partition=eap 
#SBATCH --account=Project_462000043
#SBATCH --time=10:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=4

# module load CrayEnv PrgEnv-cray craype-accel-amd-gfx90a cray-mpich rocm/5.1.4
module load cpe/22.08 PrgEnv-cray craype-accel-amd-gfx90a cray-mpich rocm/5.0.2

# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export MPICH_GPU_SUPPORT_ENABLED=1
export LD_LIBRARY_PATH=$HIP_LIB_PATH:$LD_LIBRARY_PATH


srun ./sun_realtime

