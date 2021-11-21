#!/bin/bash

#SBATCH --job-name onsim-test
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --account=Project_462000010

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128

#SBATCH -e error.out
#SBATCH -o output.out

module load cray-fftw craype-hugepages2M

touch sim_result
srun ./onsim >> sim_result

