# Makefile for phase field using p3dfft 
#
#
## # makefile for timken MEUMAPPS-SS project - martensite
FC=mpif90
FFLAGS  = -O2 -free -80
LINPLIB = $(HOME)/liblinpack.a
FFTW =  -I /nopt/nrel/apps/libraries/08-23/spack/opt/spack/linux-rhel8-icelake/intel-2021.10.0/fftw-3.3.10-aix7ofi2nlmjssdhkuaiwudzrdex2gkc/include -L /nopt/nrel/apps/libraries/08-23/spack/opt/spack/linux-rhel8-icelake/intel-2021.10.0/fftw-3.3.10-aix7ofi2nlmjssdhkuaiwudzrdex2gkc/lib -lfftw3 -lfftw3_mpi
P3DFFT  = -I /projects/ccpcmornl/for_james/P3DFFT-2.7.9/include -L /projects/ccpcmornl/for_james/P3DFFT-2.7.9/lib -lp3dfft

all: clean var_diff.x polycryst.x

OBJECTS = \
ran_2.o f_trans.o  inv_trans.o  k_space.o \
#
k_space.o : k_space.f90 
	$(FC) $(FFLAGS) -c k_space.f90  $(FFTW) $(P3DFFT) $(LINPLIB)

f_trans.o : f_trans.f90 
	$(FC) $(FFLAGS) -c f_trans.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

inv_trans.o : inv_trans.f90 
	$(FC) $(FFLAGS) -c inv_trans.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

var_diff.o : var_diff.f90 
	$(FC) $(FFLAGS) -c var_diff.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

polycryst.o : polycryst.f90 
	$(FC) $(FFLAGS) -c polycryst.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

ran_2.o : ran_2.f90 
	$(FC) $(FFLAGS) -c ran_2.f90 $(FFTW) $(P3DFFT) $(LINPLIB)

var_diff.x : var_diff.o $(OBJECTS) 
	$(FC) -o $@ $(FFLAGS) var_diff.o $(OBJECTS) $(LINPLIB) $(P3DFFT) $(FFTW)

polycryst.x : polycryst.o $(OBJECTS) 
	$(FC) -o $@ $(FFLAGS) polycryst.o $(OBJECTS) $(LINPLIB) $(P3DFFT) $(FFTW)

clean: clean
	rm *.x
