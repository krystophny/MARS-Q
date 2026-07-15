# Makefile for MARS-F on lohan6 - M. Lanctot 2008-10-01
# Taken from:
# Makefile for MARS-F on lohan4 - H. Reimerdes
# 04-21-2006 First try
#

# Macros
FFLAGS = -OPT:Olimit=0 -m64

PROG =	 marsq.x
BUILD_DIR ?= ../build
DIR_LSODE = /lsode
SRCS_LSODE = daxpy.f   dcopy.f  dewset.f  dgbsl.f  dgesl.f \
             dprepj.f  dsolsy.f  dstode.f  dvnorm.f  iumach.f  lsode.f   xerrwd.f \
             xsetun.f  dcfode.f  ddot.f   dgbfa.f   dgefa.f  dintdy.f  dscal.f \
             dsrcom.f  dumach.f  idamax.f  ixsav.f   makefile  xsetf.f


#SRCS =   $(DIR_LSODE $(SRC_LSODE))                   \

SRCS =   ./lsode/*.f    \
         dimensim.f globalm.f gijlm.f ran_modm.f  \
         feedbackm.f   gvacuumm.f  rcomdm.f \
         convolcofm.f  resmatm.f kineticm.f torquem.f \
         adaptivem.f anisotropicm.f pams.f reorbitm.f  \
         mpienv.f ToolBox.f90 marsq.f smlpck.f feedback.f printatoh.f \
         splines.f KineticExtendAnalysis.f90  \
         kinetic.f anisotropic.f adaptive.f reorbit.f \
         ntv_lib.f ntv_pre.f common_ntv.f ntv.f torque.f \
         output_t7.f fftdriver.f fftroutines.f fftoutpt.f abortrun.f

INCS =	 cdipam.inc  newrun.inc   tent.inc vtophat.inc \
         tophat.inc  comioc.inc   integc.inc    \
         compam.inc  specmat.inc \
         vtent.inc   

MODS = 

.PRECIOUS:	$(SRCS) $(INCS) $(MODS)

OBJS =	$(SRCS:.f.f90=.o)

#LIBDIR  = -L./lsode

#INCDIR = -I/usr/local/depot/openmpi-1.6-ifort-icc/include
INCDIR = 

LIBS = \
        $(LIBDIR) \
#        -llsode

#F95 = /opt/pgi/linux86-64/10.8/bin/pgf95
#F95 = ifort
F95 = mpif90 
#F95FLAGS = -mp -O2 -r8 -tp amd64e
#LDFLAGS  = -mp -O2 -r8 -tp amd64e
#F95FLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -fopenmp -static-libgfortran  
#F95FLAGS = -O2 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fopenmp -static-libgfortran  
#F95FLAGS = -fc=ifort -O2 -r8 -openmp
F95FLAGS = -O1 -r8 -openmp
#LDFLAGS =

all: $(PROG)

.PHONY: gnu ifx nvhpc clean-objects test test-runtime check
.NOTPARALLEL: gnu ifx nvhpc

gnu:
	python3 ../tools/build_with_provenance.py --profile gnu \
		--build-dir $(BUILD_DIR)

# Source the Intel oneAPI environment before invoking this target.
ifx:
	python3 ../tools/build_with_provenance.py --profile ifx \
		--build-dir $(BUILD_DIR)

# Put the NVIDIA HPC SDK compiler and HPC-X MPI wrapper on PATH first.
nvhpc:
	python3 ../tools/build_with_provenance.py --profile nvhpc \
		--build-dir $(BUILD_DIR)

$(PROG): $(OBJS)
	$(F95) $(F95FLAGS) -o $@ $(OBJS) $(LIBS) $(INCDIR)

clean:
	rm -f $(PROG) *.o *.mod *.d  \lsode\*.o \lsode\*.mod 

# The legacy dependency graph writes compiler-specific modules and objects into
# the source tree.  Remove only those intermediates before changing compiler;
# keep the historical tracked executable untouched.
clean-objects:
	rm -f *.o *.mod *.d lsode/*.o lsode/*.mod

test:
	python -B -m unittest discover -s ../tests -p 'test_*.py'

test-runtime:
	test -n "$(MARS_EXE)"
	MARS_EXE="$(MARS_EXE)" python -B -m unittest discover \
		-s ../tests -p 'test_*.py'

check: test

tar:
	gtar zcvf ../$(PROG).tar.gz $(SRCS) $(INCS) $(MODS) makefile README

.SUFFIXES: $(SUFFIXES) .f

.f90.o:
	$(F95) $(F95FLAGS) -c $<

.f.o:
	$(F95) $(F95FLAGS) -c $<
