#!/bin/bash
#SBATCH --nodes=1
##BATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
#SBATCH --partition=debug
module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel
srun -n 1 ./paraview.x > mrun.out

