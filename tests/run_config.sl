#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=phase_field_microstructure_sim
#SBATCH --output=slurm_output/sim_output.o%j 
#SBATCH --error=slurm_output/sim_error.o%j 
#SBATCH --ntasks=4
#SBATCH --account=ccpcmornl
#SBATCH --qos=high

#SBATCH --time=1:00:00
##SBATCH --partition=debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel
srun $SRUNOPTS ../src/var_diff.x ../inputs/input.txt "$1" > "$1"mrun.out
