#!/bin/bash
#SBATCH --nodes=3
#SBATCH --job-name=phase_field_microstructure_sim
#SBATCH --output=out/output.o%j
#SBATCH --error=out/error.o%j
#SBATCH --ntasks=256
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
#SBATCH --time=12:00:00
##SBATCH --partition=debug
##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel

# polycrystal
#srun $SRUNOPTS /scratch/jroger87/variable_diff/polycryst.x > out/mrun.out
srun $SRUNOPTS /scratch/jroger87/variable_diff/var_diff.x /scratch/jroger87/variable_diff/input_wetting /scratch/jroger87/phase-field-microstructure-evolution/output/variable_diff/plane_wetting_1/ > out/mrun.out
