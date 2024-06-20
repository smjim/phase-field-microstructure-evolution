#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=tests_organizer_parametric
#SBATCH --output=slurm_output/parametric_output.o%j 
#SBATCH --error=slurm_output/parametric_error.o%j 
#SBATCH --ntasks=1
#SBATCH --account=ccpcmornl
#SBATCH --qos=high

##SBATCH --partition=debug

#SBATCH --time=1:00:00

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load python/3.11.4 

# options = ['strong_test', 'weak_test', 'phi_coeff_test', 'c_coeff_test', 'd_gb_test', 'sigma1_test', 'sigma2_test', 'ap_test', 'con_0_test', 'ic_test]

# phi_coeff_test
srun python scaling_tests.py -o ../output/tmp/test_initial_conditions -i ../inputs/config.yaml --test_type phi_coeff_test --time 1:00:00
#srun python scaling_tests.py -o ../output/parametric_study/circle/phi_coeff_test -i ../inputs/config.yaml --test_type phi_coeff_test --time 1:00:00

## c_coeff_test
#srun python scaling_tests.py -o ../output/parametric_study/circle/c_coeff_test -i ../inputs/config.yaml --test_type c_coeff_test --time 1:00:00
#
## d_gb_test
#srun python scaling_tests.py -o ../output/parametric_study/circle/d_gb_test -i ../inputs/config.yaml --test_type d_gb_test --time 1:00:00
#
## sigma1_test
#srun python scaling_tests.py -o ../output/parametric_study/circle/sigma1_test -i ../inputs/config.yaml --test_type sigma1_test --time 1:00:00
#
## sigma2_test
#srun python scaling_tests.py -o ../output/parametric_study/circle/sigma2_test -i ../inputs/config.yaml --test_type sigma2_test --time 1:00:00
#
## ap_test
#srun python scaling_tests.py -o ../output/parametric_study/circle/ap_test -i ../inputs/config.yaml --test_type ap_test --time 1:00:00
#
## con_0_test
#srun python scaling_tests.py -o ../output/parametric_study/circle/con_0_test -i ../inputs/config.yaml --test_type con_0_test --time 1:00:00
#
## ic_test
#srun python scaling_tests.py -o ../output/tmp/circle_test/ic_test -i ../inputs/config.yaml --test_type ic_test --time 1:00:00
#
## con_0_test with polycrystal configuration
#srun python scaling_tests.py -o ../output/tmp/polycrystal_circle_test/con_0_test -i ../inputs/config.yaml --test_type con_0_test --time 1:00:00
