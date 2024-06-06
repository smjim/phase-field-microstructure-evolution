#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --job-name=paraview
#SBATCH --time=2:30:00
#SBATCH --ntasks=104
#SBATCH --mem-per-cpu=2G
#SBATCH --output=paraview_out/para-%j.out
#SBATCH --error=paraview_out/para-%j.err

#SBATCH --account=ccpcmornl
#SBATCH --qos=high

##SBATCH --mail-user=jroger87@vols.utk.edu
##SBATCH --mail-type=ALL

portnum=5002
paraver=5.12.1

echo "Starting Paraview headless servers on port $portnum"

paradir=$HOME/local/ParaView-5.12.1-osmesa-MPI-Linux-Python3.10-x86_64
parabin=$paradir/bin

export LD_LIBRARY_PATH=$paradir/lib/paraview-$paraver:$LD_LIBRARY_PATH
export PATH=$parabin:$PATH
export PYTHONPATH=$paradir/lib/paraview-$paraver/site-packages:$PYTHONPATH

srun $parabin/pvserver --force-offscreen-rendering --server-port=$portnum

echo "Paraview finished at $(date '+%H:%M:%S')."
