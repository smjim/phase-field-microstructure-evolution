from datetime import datetime, date
import time
import numpy as np
import csv
import os

from run_phase_sim import *
from plot_display import *

def save_output(config, sim_path, summary_file):
    duration = analyze_phase_sim_output(sim_path) 

    # Write output to file
    line = list(config.values())
    line.extend([duration, sim_path])
    with open(summary_file, 'a', newline='') as file:
        writer=csv.writer(file)
        writer.writerow(line)

if __name__ == "__main__":
    # ---------------------------
    # Step 0: Parse arguments
    # ---------------------------
    import argparse

    parser = argparse.ArgumentParser(description="Run Phase Field simulation for given parameter combination")
    parser.add_argument("-o ", "--output_dir", nargs='?', default=None, help="output directory for tests and summary file")
    parser.add_argument("--time", default="1:00:00", help="Time limit for jobs to run")

    args = parser.parse_args()
    max_runtime = args.time
    
    if args.output_dir is not None:
        output_dir = args.output_dir
    else:
        output_dir = f'/scratch/jroger87/phase-field-microstructure-evolution/output/run_{date.today().strftime("%y%m%d%H%M")}'
        print(colors.GREEN + f'Directory not specified. Using {output_dir}.' + colors.ENDC)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # ---------------------------
    # Step 1: Define parameter set
    # ---------------------------
    # Baseline configuration
    baseline_config = {
        "N_step": 10000,
        "ifreq": 1000,
        "Nx": 96,
        "Ny": 96,
        "Nz": 96,
        "t_step": 1.e-4,
        "dims": [4, 4],
        "start": 'circle',
        "grad_coeff_phi": 1.0,
        "grad_coeff_c": 1.0,
        "d_bulk": 1.0,
        "d_gb": 10.0,
        "mob_phi": 5.0,
        "sigma_1": 0.6919,
        "sigma_2": 1.02,
        "am": 0.5,
        "ap": 2.0,
        "con_0_mat": 0.3,
        "con_0_ppt": 0.7,
        "gb_force": 0.0,
        "max_runtime": '1:00:00'
    }

    # Parameter swaps
    parameter_swaps = [ # Strong scaling test
       {"dims": [1, 1], "Nx": 96, "Ny": 96, "Nz": 96},
       {"dims": [1, 2], "Nx": 96, "Ny": 96, "Nz": 96},
       {"dims": [2, 2], "Nx": 96, "Ny": 96, "Nz": 96},
       {"dims": [2, 4], "Nx": 96, "Ny": 96, "Nz": 96},
       {"dims": [4, 4], "Nx": 96, "Ny": 96, "Nz": 96},
       {"dims": [4, 8], "Nx": 96, "Ny": 96, "Nz": 96},
       {"dims": [8, 8], "Nx": 96, "Ny": 96, "Nz": 96},
    ]
    
    # Generate configuration list with parameter swaps applied to baseline
    config_list = []

    for swap in parameter_swaps:
        new_config = baseline_config.copy()
        new_config.update(swap)
        config_list.append(new_config)
    
    # Print configuration list
    print(colors.GREEN + 'Configs to be tested:' + colors.ENDC)
    print(colors.RED)
    print_configs(config_list)
    print(colors.ENDC)

    # Write output summary header
    summary_file = os.path.join(output_dir, 'output.csv')
    with open(summary_file, 'w', newline='') as file:
        writer=csv.writer(file) 
        header=list(baseline_config.keys())         # Input parameters
        header.extend(['duration', 'output_dir'])   # Tested outputs
        writer.writerow(header)

    # ---------------------------
    # Step 2: Run simulation with parameter configurations 
    # ---------------------------
    job_dict = {}

    # Save submitted jobs to job id dict for analysis upon completion
    for i, config in enumerate(config_list):
        # Run simulation, save duration data
        sim_path, jobid = run_phase_sim(output_dir, config)
        if jobid is not None:
            job_dict[jobid] = {'config': config, 'sim_path':sim_path}
        else:   # simulation has already been run, analyze output
            save_output(config, sim_path, summary_file)

    # ---------------------------
    # Step 3: Analyze simulation outputs
    # ---------------------------
    # Track which jobs are completed
    while job_dict:
        completed_jobs = []
        print(f'Checking job status {datetime.now().strftime("%H:%M:%S")}')
        for job_id, job_info in job_dict.items():
            if check_job_status(job_id):
                # job is complete, analyze the output
                save_output(job_info['config'], job_info['sim_path'], summary_file)

        # Remove completed jobs from the dictionary
        for job_id in completed_jobs:
            del job_dict[job_id]

        # Check again after some break
        time.sleep(60)
