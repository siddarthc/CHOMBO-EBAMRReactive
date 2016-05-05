# -*- Mode: Makefile -*-

# trace the chain of include makefiles
makefiles += NASA4_exec

# This is an embedded boundary example
USE_EB=TRUE

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = ebamrReactive

# the location of the Chombo "lib" directory
#CHOMBO_HOME = ../../lib
CHOMBO_HOME = /work/02138/siddarth/CHOMBO/EBAMRReactive/lib


# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRTimeDependent EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTimeDependent AMRTools BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../src ../chemkin

# input file for 'run' target
INPUT = ramp.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets



