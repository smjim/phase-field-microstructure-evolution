#!/bin/bash
#SBATCH --job-name=paraview_visualization
#SBATCH --output=slurm_output/output.o%j 
#SBATCH --error=slurm_output/error.o%j 
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1:00:00
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
##SBATCH -p debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

# Create .vtk data
module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel

# DIR="out/"
# mkdir "${DIR}visualization"
# srun -n 4 /scratch/jroger87/phase-field-microstructure-evolution/src/paraview.x "$DIR" "$DIR"visualization 0.001 > "$DIR"visualization/precipitate_fraction.dat

# DIR="/scratch/jroger87/phase-field-microstructure-evolution/output/variable_diff/plane_output_nowetting/"
# OUTDIR="/scratch/jroger87/phase-field-microstructure-evolution/output/variable_diff/plane_output_nowetting/output/"

DIR="/scratch/jroger87/phase-field-microstructure-evolution/output/variable_diff/plane_wetting_1/"
OUTDIR="/scratch/jroger87/phase-field-microstructure-evolution/output/variable_diff/plane_wetting_1/visualization/"

srun -n 16 /scratch/jroger87/phase-field-microstructure-evolution/src/paraview.x "$DIR" "$OUTDIR" 0.001 > "$OUTDIR"precipitate_fraction.dat
