#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=parallel_weak_tests_organizer
#SBATCH --output=slurm_output/parallel_weak_output.o%j 
#SBATCH --error=slurm_output/parallel_weak_error.o%j 
#SBATCH --ntasks=1
#SBATCH --account=ccpcmornl
#SBATCH --qos=high

#SBATCH --time=5:00:00
##SBATCH --partition=debug

#SBATCH --mail-user=jroger87@vols.utk.edu
#SBATCH --mail-type=ALL

module purge
module load python/3.11.4 

# Iteration 0
srun python scaling_tests.py -o ../output/parallel_weak/iteration_0/ -i ../inputs/config.yaml --test_type weak_test --time 0:45:00 
## Iteration 1
#srun python scaling_tests.py -o ../output/parallel_weak/iteration_1/ -i ../inputs/config.yaml --test_type weak_test --time 0:45:00 
# Iteration 2
srun python scaling_tests.py -o ../output/parallel_weak/iteration_2/ -i ../inputs/config.yaml --test_type weak_test --time 0:45:00 
# Iteration 3
srun python scaling_tests.py -o ../output/parallel_weak/iteration_3/ -i ../inputs/config.yaml --test_type weak_test --time 0:45:00 
# Iteration 4
srun python scaling_tests.py -o ../output/parallel_weak/iteration_4/ -i ../inputs/config.yaml --test_type weak_test --time 0:45:00 
# Iteration 5
srun python scaling_tests.py -o ../output/parallel_weak/iteration_5/ -i ../inputs/config.yaml --test_type weak_test --time 0:45:00 
