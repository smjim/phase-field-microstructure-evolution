#!/bin/bash
#SBATCH --job-name=plane_analysis
#SBATCH --output=slurm_output/output.o%j 
#SBATCH --error=slurm_output/error.o%j 
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=0:10:00
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
##SBATCH -p debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

# Create .vtk data
module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel

DIR="/scratch/jroger87/variable_diff/"
OUTDIR="/scratch/jroger87/phase-field-microstructure-evolution/output/variable_diff/plane_wetting_1/"
srun -n 4 /scratch/jroger87/phase-field-microstructure-evolution/src/analyze_plane.x "$DIR" "$OUTDIR" 0.1 
