#!/bin/bash
#SBATCH --job-name=paraview_visualization
#SBATCH --output=slurm_output/output.o%j 
#SBATCH --error=slurm_output/error.o%j 
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=0:15:00
#SBATCH --account=ccpcmornl
#SBATCH --qos=high
##SBATCH -p debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

# Create .vtk data
module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel

TMPDIR=""
srun -n 4 ../src/paraview.x "$TMPDIR" "$TMPDIR"visualization 0.2 > "$TMPDIR"visualization/precipitate_fraction.dat

# process all "runs" within all "outdirs" within a containing directory
ININDIR="../output/wetting_study/polycrystal/"
for INDIR in "$ININDIR"*/; do
    for DIR in "$INDIR"*/; do
        if [[ -d "$DIR" ]]; then
            echo "Processing directory: $DIR"
            echo "srun -n 4 ../src/paraview.x $DIR $DIR visualization 0.2"
            srun -n 4 ../src/paraview.x "$DIR" "$DIR"visualization 0.2 > "$DIR"visualization/precipitate_fraction.dat
        fi
    done
done

# # process all "runs" within an "outdir" 
# INDIR="/scratch/jroger87/phase-field-microstructure-evolution/output/wetting_study/polycrystal/d_gb_test/"
# for DIR in "$INDIR"*/; do
#     if [[ -d "$DIR" ]]; then
#         echo "Processing directory: $DIR"
#         echo "srun -n 4 ../src/paraview.x $DIR $DIR visualization 0.2"
#         mkdir "$DIR"visualization
#         srun -n 4 ../src/paraview.x "$DIR" "$DIR"visualization 0.2 > "$DIR"visualization/precipitate_fraction.dat
#     fi
# done
