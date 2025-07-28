# Makefile for MARS-F on lohan4 - H. Reimerdes
# 04-21-2006 First try
#

# Macros
FFLAGS = -OPT:Olimit=0 -m64

PROG =	 chease.x

SRCS =   chease.f 

INCS =	       COMBND.inc  COMINT.inc  COMPLO.inc  DECLAR.inc	QUAQDQ.inc \
               COMCON.inc  COMIOD.inc  COMRCH.inc  HERMIT.inc	QUAQQQ.inc \
               BNDIND.inc  COMDAT.inc  COMISO.inc  COMSOL.inc           \
               COMDIM.inc  COMLAB.inc  COMSUR.inc  SOLOV.inc            \
               COMEQD.inc  COMMAP.inc  COMVAC.inc  NEWRUN.inc	        \
               COMERA.inc  COMNUM.inc  COMVEV.inc                       \
               COMBAL.inc  COMESH.inc  COMOPT.inc  CUCCCC.inc           \
               COMBLA.inc     COMETA.inc  COMPHY.inc  CUCDCD.inc  PROCESS.inc

MODS = 

.PRECIOUS:	$(SRCS) $(INCS) $(MODS)

OBJS =	$(SRCS:.f=.o)

#LIBDIR = \
#        -L/usr/pgi/linux86-64/7.0-5/lib

#LIBS = \
#        $(LIBDIR)
LIBS =

#F95 = pgf95
#F95 = ifort 
F95 = mpif90
 
#F95FLAGS = -O2 -r8 -tp amd64e
#F95FLAGS = -O2 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fopenmp -static-libgfortran  
F95FLAGS = -O3 -r8 
LDFLAGS =

all: $(PROG)

$(PROG): $(OBJS)
	$(F95) $(F95FLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.d work.*

test:
	echo $(OBJS) 

tar:
	gtar zcvf ../$(PROG).tar.gz $(SRCS) $(INCS) $(MODS) Makefile_Portland README

.SUFFIXES: $(SUFFIXES) .f

.f.o:
	$(F95) $(F95FLAGS) -c $<

