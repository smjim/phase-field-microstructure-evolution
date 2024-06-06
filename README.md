# phase-field-microstructure-evolution
Phase field simulation of microstructure evolution

## Table of Contents
1. Description   
	1.1 Theory   
	1.2 Features    
2. Installation  
3. Usage  
4. Author(s) and Acknowledgement  
5. References  

## 1 Description

### 1.1 Theory of Phase-Field Microstructure Evolution Simulation 

### 1.2 Features

## 2 Installation

## 3 Usage
1. Compile `var_diff.x` and `paraview.x` using `$ make` in`src/`
2. Determine parameter default settings and swaps by modifying configurations in `inputs/config.yaml`
3. Configure slurm scripts and computer/ user specific parameters in `tests/run_phase_sim.py` and `tests/scaling_tests.py`
4. Run `$ python scaling_tests.py -o <run_directory> -i <input_config_yaml> --test_type <test_type> --time <max_time_per_run>` from within `tests/`  
 - Make sure to add the option `--visualize` if you wish for `paraview.x` to automatically create .vtk files from the simulation output
5. Visualize results with Paraview after conversion to .vtk file format  
 - This can be done with the `paraview_headless.sl` script, after installing Paraview headless on server and Paraview client on PC, and by connecting to server via specified port (for example by piping port 5002 traffic to port 5002 on local computer).

## 4 Authors and Acknowledgement
James M. Rogers, University of Tennessee - Knoxville  
contact: jroger87@vols.utk.edu  

## 5 References

## License  

## Project Status  
Currently in progress as of 6-6-2024 
