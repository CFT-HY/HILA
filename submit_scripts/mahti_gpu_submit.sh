#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --partition=gpumedium
#SBATCH --account=Project_2005155

#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32

#SBATCH -J onsim_large
#SBATCH --gres=gpu:a100:4
# ## #SBATCH --gpu-bind=single:1

### #SBATCH --hint=nomultithread
#SBATCH -e error.out
#SBATCH -o out

ml purge
ml gcc/10.3.0 cuda/11.4.2 openmpi/4.1.0

touch sim_result
UCX_MEMTYPE_CACHE=n HCOLL_GPU_CUDA_MEMTYPE_CACHE_ENABLE=0 HCOLL_GPU_ENABLE=1 srun ./onsim >> sim_result

