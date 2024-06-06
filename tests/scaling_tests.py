from datetime import date
import yaml
import time
import numpy as np
import csv
import os

from run_phase_sim import *
from plot_display import *

def load_config(yaml_file):
    with open(yaml_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

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
    parser.add_argument("-i ", "--input_params", nargs='?', default='../inputs/config.yaml', help="input .yaml file for parameter combinations to be tested")
    parser.add_argument("-o ", "--output_dir", nargs='?', default=None, help="output directory for tests and summary file")
    parser.add_argument("--test_type", type=str, choices=['strong_test', 'weak_test', 'phi_coeff_test', 'c_coeff_test', 'd_gb_test', 'sigma1_test', 'sigma2_test', 'ap_test', 'con_0_test', 'ic_test'], required=False, help="Time limit for jobs to run")
    parser.add_argument("--time", default="1:00:00", help="Time limit for jobs to run")
    parser.add_argument('--visualize', action='store_true', help='Automatically generate paraview visualization .vtk files')

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
    config = load_config(args.input_params)
    baseline_config = config['baseline_config']

    # Choose parameter set
    if args.test_type is not None:
        parameter_swaps = config[args.test_type]    # Specify via command line
    else:
        parameter_swaps = config['weak_test']       # Manually specify
    
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
        sim_path, jobid = run_phase_sim(output_dir, config, max_runtime)
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

        # Boolean array of completion status
        job_status = check_job_status(list(job_dict.keys()))

        for job_id, job_info in job_dict.items():
            if job_id in job_status:
                # job is complete, analyze the output
                save_output(job_info['config'], job_info['sim_path'], summary_file)

                # run the visualization creation tool
                if args.visualize:
                    run_paraview_converter(job_info['sim_path'])

                completed_jobs.append(job_id)  # Add completed job to the list        

        # Remove completed jobs from the dictionary
        for job_id in completed_jobs:
            del job_dict[job_id]

        # Check again after some break
        time.sleep(30)
