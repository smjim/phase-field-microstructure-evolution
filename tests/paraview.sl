#!/bin/bash
#SBATCH --job-name=paraview_visualization
#SBATCH --output="$1"output.o%j 
#SBATCH --error="$1"error.o%j 
#SBATCH --nodes=1
##BATCH --ntasks=1
#SBATCH --time=0:05:00
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
#SBATCH -p debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

# Create .vtk data
module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel
srun -n 1 ../src/paraview.x "$1" "$1"visualization/ > "$1"visualization/mrun.out
