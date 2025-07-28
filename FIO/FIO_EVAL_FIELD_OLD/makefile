# Macros
FFLAGS = -OPT:Olimit=0 -m64

PROG =	 fio_test.x

SRCS =   fio_eval_field_fm.f fio_eval_field_f.f fio_test.f 

.PRECIOUS:	$(SRCS) 

OBJS =	$(SRCS:.f.f90=.o)

LIBS = \
        $(LIBDIR) \

F95      = mpif90 
F95FLAGS = -O1 -r8 -openmp

all: $(PROG)

$(PROG): $(OBJS)
	$(F95) $(F95FLAGS) -o $@ $(OBJS) $(LIBS) 

clean:
	rm -f $(PROG) *.o *.mod *.d  

test:
	echo $(OBJS) 

.SUFFIXES: $(SUFFIXES) .f

.f90.o:
	$(F95) $(F95FLAGS) -c $<

.f.o:
	$(F95) $(F95FLAGS) -c $<


