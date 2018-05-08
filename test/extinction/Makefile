VPATH=../../src
obj=exponential_integral.o extinction.o test.o

default: $(obj)
	gfortran -g -O0 --backtrace $(obj) -o test

%.o:%.f90
	gfortran -g -O0 --backtrace -c $^ -o $@

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec)