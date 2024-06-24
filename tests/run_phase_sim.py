from datetime import datetime
import numpy as np
import subprocess
import hashlib
import csv 
import re
import os

# Define the base SLURM script template
slurm_template = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=phase_field_microstructure_sim
#SBATCH --output={output_path}
#SBATCH --error={error_path}
#SBATCH --ntasks={ntasks}
#SBATCH --account=ccpcmornl
##SBATCH --qos=high
#SBATCH --time={time}
##SBATCH --partition=debug
##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel

srun $SRUNOPTS /scratch/jroger87/phase-field-microstructure-evolution/src/var_diff.x {input_file_path} {output_dir} > {mrun_path}
"""

paraview_slurm_template = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=paraview_converter
#SBATCH --output={output_file}
#SBATCH --error={error_file}
#SBATCH --ntasks={ntasks}
##SBATCH --qos=high
#SBATCH --account=ccpcmornl
#SBATCH --time={time}
##SBATCH --partition=debug
##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

module purge
module load PrgEnv-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel

srun $SRUNOPTS {fortran_executable} {sim_dirname} {vis_outdir} {threshold_str} > {precipitate_outfile}
"""

class colors:
    RED = '\033[31m'
    ENDC = '\033[m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'

# for unique naming of output files
def generate_dirname(output_dir, parameters):
    # Define a regular expression pattern to match allowed characters (alphanumeric and underscores)
    pattern = re.compile(r'[^a-zA-Z0-9_\.]')

    # Process parameters and replace unwanted characters
    parameterlist = []
    for value in parameters.values():
        if isinstance(value, list):
            value_str = '_'.join(map(str, value))
        else:
            value_str = str(value)
        # Replace unwanted characters with underscores
        value_str = pattern.sub('_', value_str)
        parameterlist.append(value_str)
    
    combined_parameters = '_'.join(parameterlist)

    ## Parameter dirname
    #output_dirname = os.path.join(output_dir, f'run_{combined_parameters}')

    # Hash dirname
    hash_object = hashlib.sha256(combined_parameters.encode())
    hash_value = hash_object.hexdigest()
    output_dirname = os.path.join(output_dir, f'run_{hash_value}')

    return output_dirname

def check_job_status(job_ids):
    print(f'Checking job status {datetime.now().strftime("%H:%M:%S")}')

    try:
        # Print the output status of all job_ids
        full_result = subprocess.run(["squeue", "-u", "jroger87"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        print(full_result.stdout)

        # Run squeue with options and get the output
        result = subprocess.run(["squeue", "--noheader", "-o", "%i"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        
        #print("Result stdout:", result.stdout)
        #print("Result stderr:", result.stderr)
    
        # Extract job IDs from the squeue output
        running_job_ids = set(result.stdout.split())
    
        # Determine which jobs are complete (not running or pending)
        completed_jobs = [job_id for job_id in job_ids if job_id not in running_job_ids]

        return completed_jobs

    except subprocess.CalledProcessError as e:
        print(f"Error checking job status: {e.stderr}")
        return []

# generate input file
def generate_input_file(output_dir, params):
    input_file_path = os.path.join(output_dir, 'input.txt')

    with open(input_file_path, 'w') as file:
        file.write(f"1       / nrun\n")
        file.write(f"4 4 4 {params['N_step']} {params['ifreq']}      / ppt_rad(3), N_step, ifreq\n")
        file.write(f"48 70 48 -1        / i_ppt, j_ppt, k_ppt, iseed\n")
        file.write(f"40 {params['num_ppt']}  / grn_rad, num_ppt\n")
        file.write(f"{params['Nx']} {params['Ny']} {params['Nz']} {params['var']} 0.01 0.01 0.01 {params['t_step']}  / Nx, Ny, Nz, var, dx, dy, dz, t_step\n")
        file.write(f"{params['grad_coeff_phi']} {params['grad_coeff_c']} {params['d_bulk']} {params['d_gb']} {params['mob_phi']}   / grad_coeff_phi, grad_coeff_c, d_bulk, d_gb, mob_phi\n")
        file.write(f"6.0  3.0  / sigma_1, sigma_2\n") 
        file.write(f"{params['am']} {params['ap']} {params['con_0_mat']} {params['con_0_ppt']} {params['gb_force']}     / am, ap, con_0_mat, con_0_ppt, gb_force\n")
        file.write(f"2       / ndim\n")
        file.write(f"{params['dims'][0]} {params['dims'][1]}     / dims(2)\n")
        file.write(f"'{params['start']}'         / start\n")
        file.write(f"{params['sigma_1']}  {params['sigma_2']}  / sigma_1, sigma_2\n")
    
    return input_file_path
    
def submit_slurm_job(config, output_dir, input_file_path, time_str="1:00:00"):
    # Calculate the number of tasks
    ntasks = int(config["dims"][0] * config["dims"][1])

    # Determine additional file paths
    mrun_path = os.path.join(output_dir, 'mrun.out')
    output_path = os.path.join(output_dir, 'output.o%j')
    error_path = os.path.join(output_dir, 'error.o%j')
    
    # Create the SLURM script content
    slurm_content = slurm_template.format(
        output_dir=f'{output_dir}/',
        ntasks=ntasks,
        time=time_str,
        mrun_path=mrun_path,
        output_path=output_path,
        input_file_path=input_file_path,
        error_path=error_path
    )
    
    # Define the SLURM script file name
    script_filename = os.path.join(output_dir, "slurm_job.sh")
    
    # Write the SLURM script to a file
    with open(script_filename, 'w') as script_file:
        script_file.write(slurm_content)
    
    # Submit the SLURM job using sbatch
    result = subprocess.run(["sbatch", script_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    
    print("Result stdout:", result.stdout)
    print("Result stderr:", result.stderr)
    
    if result.returncode != 0:
        print(f"Error submitting SLURM job: {result.stderr}")
        return None, None

    # Extract job ID from the command output
    if result.stdout:
        job_id = result.stdout.strip().split()[-1]
        print(f"SLURM job submitted successfully. Job ID: {job_id}")
    else:
        print("No output captured from sbatch command.")
        job_id = None

    return script_filename, job_id

# run simulation with given parameter combination
def run_phase_sim(output_dir, params, max_runtime):
    # If simulation has been run previously, dont rerun it
    sim_dirname = generate_dirname(output_dir, params)
    if not os.path.exists(sim_dirname):
        os.makedirs(sim_dirname)
    else:
        print(f'Simulation has already been run: {sim_dirname}')
        return sim_dirname, None

    # Generate input.txt file according to parameters
    input_file_path = generate_input_file(sim_dirname, params) 
    print(colors.GREEN + "Using input file path:" + colors.ENDC)
    print(colors.BLUE + input_file_path + colors.ENDC)
    
    # Run var_diff.x with slurm according to parameters
    slurm_filename, jobid = submit_slurm_job(params, sim_dirname, input_file_path, max_runtime) # run the job
    #slurm_filename, jobid = os.path.join(sim_dirname, 'tmp.sl'), None                          # dont actually run the job
    print(colors.GREEN + "Running Slurm File:" + colors.ENDC)
    print(colors.BLUE + slurm_filename + colors.ENDC)

    return sim_dirname, jobid
    
# analyze simulation output for duration 
def analyze_phase_sim_output(sim_dirname):
    # Calculate computation duration of run
    def calculate_duration(sim_dirname):
        # Duration comes from os.path.join(sim_dirname, 'time_step.dat')
        duration = None
        with open(os.path.join(sim_dirname, 'time_step.dat'), 'r') as file:
            reader = csv.reader(file, delimiter=' ', skipinitialspace=True)
            next(reader)  # Skip header line
            for row in reader:
                if len(row) == 0:  # Skip empty lines
                    continue
                duration = float(row[2])
        return duration

    # Calculate fraction of boundaries covered by precipitates with fortran
    def calculate_final_phase_fraction(sim_dirname):
        threshold = 0.2
        # Run the visualization creation tool
        phase_fraction_file = run_paraview_converter(sim_dirname, threshold)
    
        # Read the last line from the output file
        if phase_fraction_file is not None:
            with open(phase_fraction_file, 'r') as file:
                lines = file.readlines()
                if not lines:
                    raise ValueError(f"Output file {phase_fraction_file} is empty.")
                last_line = lines[-1].strip()
                phase_fraction = float(last_line.split()[3])  # Assuming the phase fraction is in the fourth column

            return phase_fraction
        else:
            return "rerun for analysis"

    # ------------------------------------
    # TODO add other quantities to be measured
    # ------------------------------------

    # Calculate duration of simulation runtime
    duration = calculate_duration(sim_dirname)

    # Calculate precipitate coverage of boundaries
    phase_fraction = None #calculate_final_phase_fraction(sim_dirname)

    # ------------------------------------

    return duration, phase_fraction

# Generate .vtk files from step_ files to visualize result with paraview
def run_paraview_converter(sim_dirname, threshold):

    # Create visualization directory if not exists
    vis_outdir = os.path.join(sim_dirname, 'visualization')

    if not os.path.exists(vis_outdir):
        os.makedirs(vis_outdir)

    # Construct the command to run the Fortran program
    fortran_executable = "../src/paraview.x"
    precipitate_outfile = os.path.join(vis_outdir, 'precipitate_fraction.dat')
    output_file = os.path.join(vis_outdir,'output.o%j')
    error_file = os.path.join(vis_outdir,'error.o%j')

    #vis_outdir=f'{vis_outdir}/',
    slurm_content = paraview_slurm_template.format(
        vis_outdir=vis_outdir,
        sim_dirname=sim_dirname,
        ntasks=4,
        time='0:10:00',
        fortran_executable=fortran_executable,
        precipitate_outfile=precipitate_outfile,
        output_file = output_file,
        error_file = error_file,
        threshold_str=str(threshold)
    )
    
    # Define the SLURM script file name
    script_filename = os.path.join(vis_outdir, "slurm_job.sh")
    command = f'srun --time=0:10:00 --account=ccpcmornl -n 4 {fortran_executable} {sim_dirname} {vis_outdir} {str(threshold)} > {precipitate_outfile}'
    print(colors.GREEN + "Running slurm file:" + colors.ENDC)
    print(colors.BLUE + script_filename + colors.ENDC)
    print(colors.GREEN + "Containing command:" + colors.ENDC)
    print(colors.BLUE + command + colors.ENDC)
    
    # Write the SLURM script to a file
    with open(script_filename, 'w') as script_file:
        script_file.write(slurm_content)
    
    # Submit the SLURM job using sbatch
    result = subprocess.run(["sbatch", script_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    print(colors.GREEN + "Slurm contents:" + colors.ENDC)
    print(colors.RED + slurm_content + colors.ENDC)
    
    if result.returncode != 0:
        print("Result stdout:", result.stdout)
        print("Result stderr:", result.stderr)
    
        print(f"Error submitting SLURM job: {result.stderr}")
        print(colors.GREEN + "Slurm contents:" + colors.ENDC)
        print(colors.RED + slurm_content + colors.ENDC)
        return None
    
    # Check if the output file exists
    if not os.path.exists(precipitate_outfile):
        return None 
        #raise FileNotFoundError(f"Output file {precipitate_outfile} not found.")

    return precipitate_outfile 
