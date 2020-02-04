# Hans A. Winther (2015) (hans.a.winther@gmail.com)

SHELL := /bin/bash

# Set compiler
CC = g++ -std=c++11 

# Paths to GSL library
INC  = -I/Users/jakobs/c++_stuff/gsl-2.6/include
#INC += -I/Users/jakobs/miniconda3/include/python3.7m
#INC += -I/usr/local/include/python2.7
LIBS = -L/Users/jakobs/c++_stuff/gsl-2.6/lib -lgsl -lgslcblas
#LIBS += -L/Users/jakobs/miniconda3/include/python3.7m# -lpython3.7m

#=======================================================
# Options
#=======================================================
OPTIONS = 

# Add bounds checking
OPTIONS += -D_GLIBCXX_DEBUG

# Show warnings if atempting to evaluate a spline out of bounds
OPTIONS += -D_SPLINE_WARNINGS_ON

# Show info about the solution as we integrate
# OPTIONS = -D_FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE

# Add OpenMP parallelization
# OPTIONS += -D_USEOPEMP
# CC += -fopenmp

# Add bessel function library (otherwise use the GSL one)
#OPTIONS += -D_COMPLEX_BESSEL
#LIBS += -lgfortran -lcomplex_bessel
#=======================================================

C = -O3 $(OPTIONS)

TARGETS := cmb
all: $(TARGETS)

# OBJECT FILES
OBJS = Main.o Utils.o BackgroundCosmology.o RecombinationHistory.o Perturbations.o PowerSpectrum.o Spline.o ODESolver.o writeme.o

# DEPENDENCIES
Main.o		            : BackgroundCosmology.h RecombinationHistory.h Perturbations.h PowerSpectrum.h
Spline.o                : Spline.h
ODESolver.o             : ODESolver.h
Utils.o                 : Utils.h Spline.h ODESolver.h
BackgroundCosmology.o	: BackgroundCosmology.h Utils.h Spline.h ODESolver.h
RecombinationHistory.o  : RecombinationHistory.h BackgroundCosmology.h
Perturbations.o         : Perturbations.h BackgroundCosmology.h RecombinationHistory.h
PowerSpectrum.o         : PowerSpectrum.h BackgroundCosmology.h RecombinationHistory.h Perturbations.h
Examples.o              : Utils.h Spline.h ODESolver.h writeme.h
writeme.o  				: writeme.h

examples: Examples.o Utils.o Spline.o ODESolver.o writeme.o
	${CC} -o $@ $^ $C $(INC) $(LIBS)

cmb: $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIBS)

%.o: %.cpp
	${CC}  -c -o $@ $< $C $(INC) 

clean:
	rm -rf $(TARGETS) *.o

