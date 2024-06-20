#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=generate_vtk
#SBATCH --output=slurm_output/genvtk_output.o%j 
#SBATCH --error=slurm_output/genvtk_error.o%j 
#SBATCH --ntasks=1
#SBATCH --account=ccpcmornl
#SBATCH --qos=high

#SBATCH --time=1:00:00
#SBATCH --partition=debug

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load python/3.11.4 

# Specify input directory 
srun python generate_vtk.py -i ../output/parametric_study/polycrystal/ic_test/
