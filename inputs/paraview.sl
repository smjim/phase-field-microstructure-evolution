#!/bin/bash
#SBATCH --job-name=paraview_visualization
#SBATCH --output=output/output.o%j 
#SBATCH --error=error/error.o%j 
#SBATCH --nodes=1
##BATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
#SBATCH --partition=debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel
srun -n 1 ./paraview.x $1 > ../outputs/mrun.out
