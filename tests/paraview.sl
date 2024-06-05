#!/bin/bash
#SBATCH --job-name=paraview_visualization
#SBATCH --output=../output/run_5-31/visualization/output.o%j 
#SBATCH --error=../output/run_5-31/visualization/error.o%j 
#SBATCH --nodes=1
##BATCH --ntasks=1
#SBATCH --time=0:05:00
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
#SBATCH --partition=debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel
srun -n 1 ../src/paraview.x ../output/run_5-31/ ../output/run_5-31/visualization/ > ../output/run_5-31/visualization/mrun.out
