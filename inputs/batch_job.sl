#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=phase_field_microstructure_sim
#SBATCH --output=../output/output.o%j 
#SBATCH --error=../error/error.o%j 
#SBATCH --ntasks=16
#SBATCH --account=ccpcmornl
#SBATCH --qos=high

#SBATCH --time=2:00:00
##SBATCH --partition=debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel
srun $SRUNOPTS ../src/var_diff.x ../output/ > ../output/mrun.out
