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
module load python/3.11.4

srun python scaling_tests.py -o ../output/tmp/polycrystal_test -i ../inputs/config.yaml --test_type phi_coeff_test --time 1:00:00
