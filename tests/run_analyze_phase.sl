#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=analyze_vtk
#SBATCH --output=slurm_output/vtk_output.o%j 
#SBATCH --error=slurm_output/vtk_error.o%j 
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

# `$ ./analyze_phase.x <input_dir>`
#srun $SRUNOPTS ../src/analyze_phase.x ../output/parametric_study/plane/phi_coeff_test/run_8a1f862051b45079b364585780f443d0e35b915a95b69dc5de214acf8a3801e9 ../output/parametric_study/plane/phi_coeff_test/run_8a1f862051b45079b364585780f443d0e35b915a95b69dc5de214acf8a3801e9/tmp
srun $SRUNOPTS ../src/analyze_phase.x ../output/parametric_study/circle/phi_coeff_test/run_c4d97e424dcc02b44d9e77806b2a91133eebad7b80bd53425a639572b6aeb3e4 ../output/parametric_study/circle/phi_coeff_test/run_c4d97e424dcc02b44d9e77806b2a91133eebad7b80bd53425a639572b6aeb3e4/visualization 

