#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=phase_field_microstructure_sim
#SBATCH --output=out1/output.o%j
#SBATCH --error=out1/error.o%j
#SBATCH --ntasks=64
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
#SBATCH --time=2:00:00
##SBATCH --partition=debug
##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel

# circle
#srun $SRUNOPTS /scratch/jroger87/phase-field-microstructure-evolution/src/var_diff.x /scratch/jroger87/phase-field-microstructure-evolution/src/tmp/input_circle out/ > out/mrun.out

# polycrystal bonus
srun $SRUNOPTS /scratch/jroger87/phase-field-microstructure-evolution/src/var_diff.x /scratch/jroger87/phase-field-microstructure-evolution/src/tmp/input out/ > out/mrun.out
#srun $SRUNOPTS /scratch/jroger87/phase-field-microstructure-evolution/src/var_diff.x /scratch/jroger87/phase-field-microstructure-evolution/src/tmp/input_poly_bonus out/ > out/mrun.out
#srun $SRUNOPTS /scratch/jroger87/phase-field-microstructure-evolution/src/var_diff.x /scratch/jroger87/phase-field-microstructure-evolution/src/tmp/input_polycrystal out/ > out/mrun.out

## polycrystal
#srun $SRUNOPTS /scratch/jroger87/phase-field-microstructure-evolution/src/var_diff.x /scratch/jroger87/phase-field-microstructure-evolution/src/tmp/input_polycrystal out/ > out/mrun.out

## polycrystal with valgrind
#srun $SRUNOPTS valgrind --leak-check=yes --track-origins=yes /scratch/jroger87/phase-field-microstructure-evolution/src/var_diff.x /scratch/jroger87/phase-field-microstructure-evolution/src/tmp/input_polycrystal out/ > out/mrun.out

# plane bonus
#srun $SRUNOPTS /scratch/jroger87/phase-field-microstructure-evolution/src/var_diff.x /scratch/jroger87/phase-field-microstructure-evolution/src/tmp/input_plane out/ > out/mrun.out
