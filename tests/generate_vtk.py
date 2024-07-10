import argparse
import os

from run_phase_sim import run_paraview_converter

def generate_vtk_files(input_dir):
    # List all directories in the input directory
    for subdir in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, subdir)
        # Check if the subdirectory name matches the pattern 'run_*' and is a directory
        if os.path.isdir(subdir_path) and subdir.startswith('run_'):
            print(f"Processing directory: {subdir_path}")
            run_paraview_converter(subdir_path)

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="Generate .vtk files for given input directory")
    parser.add_argument("-i", "--input_dir", required=False, help="Input directory for run directories")
    parser.add_argument("--run", required=False, help="Run directory for direct paraview conversion")

    args = parser.parse_args()

    if args.input_dir:
        # Generate .vtk files for the directories in the input directory
        generate_vtk_files(args.input_dir)

    elif args.run:
        run_paraview_converter(args.run)
