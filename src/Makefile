# Makefile for phase field using p3dfft 
# timken MEUMAPPS-SS project - martensite

FC=mpif90
FFLAGS  = -O2 -free -80 # Original
#FFLAGS  = -O0 -free -80 -g # Modified for debugging
LINPLIB = $(HOME)/liblinpack.a
FFTW =  -I /nopt/nrel/apps/libraries/08-23/spack/opt/spack/linux-rhel8-icelake/intel-2021.10.0/fftw-3.3.10-aix7ofi2nlmjssdhkuaiwudzrdex2gkc/include -L /nopt/nrel/apps/libraries/08-23/spack/opt/spack/linux-rhel8-icelake/intel-2021.10.0/fftw-3.3.10-aix7ofi2nlmjssdhkuaiwudzrdex2gkc/lib -lfftw3 -lfftw3_mpi
P3DFFT  = -I /projects/ccpcmornl/for_james/P3DFFT-2.7.9/include -L /projects/ccpcmornl/for_james/P3DFFT-2.7.9/lib -lp3dfft

all: clean var_diff.x paraview.x analyze_plane.x

OBJECTS = \
  ran_2.o f_trans.o inv_trans.o k_space.o #initialize_conditions.o

#MODULES = initialize_conditions.o

var_diff.x: var_diff.o $(OBJECTS)
	$(FC) -o $@ $(FFLAGS) var_diff.o $(OBJECTS) $(LINPLIB) $(P3DFFT) $(FFTW)

k_space.o: k_space.f90
	$(FC) $(FFLAGS) -c k_space.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

f_trans.o: f_trans.f90
	$(FC) $(FFLAGS) -c f_trans.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

inv_trans.o: inv_trans.f90
	$(FC) $(FFLAGS) -c inv_trans.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

var_diff.o: var_diff.f90 $(MODULES)
	$(FC) $(FFLAGS) -c var_diff.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

ran_2.o: ran_2.f90
	$(FC) $(FFLAGS) -c ran_2.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

#initialize_conditions.o: initialize_conditions.f90 ran_2.o
#	$(FC) $(FFLAGS) -c initialize_conditions.f90 ran_2.o $(FFTW) $(P3DFFT) $(LINPLIB)

# Paraview visualization and analysis of phase data
paraview.x: paraview.o ran_2.o 
	$(FC) -o $@ $(FFLAGS) paraview.o ran_2.o 

paraview.o: paraview.f90
	$(FC) $(FFLAGS) -c paraview.f90

# Analysis of plane wetting
analyze_plane.x: analyze_plane.o ran_2.o
	$(FC) -o $@ $(FFLAGS) analyze_plane.o ran_2.o

analyze_plane.o: analyze_plane.f90
	$(FC) $(FFLAGS) -c analyze_plane.f90

# Analysis of polycryst wetting
analyze_polycryst.x: analyze_polycryst.o ran_2.o
	$(FC) -o $@ $(FFLAGS) analyze_polycryst.o ran_2.o

analyze_polycryst.o: analyze_polycryst.f90
	$(FC) $(FFLAGS) -c analyze_polycryst.f90

clean:
	rm -f *.o *.mod var_diff.x paraview.x analyze_plane.x
