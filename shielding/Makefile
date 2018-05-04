# Makefile for test_leiden_shielding.f90

fc=gfortran

switch=-ffree-line-length-none
nowarn=-w

#objects
objs = opkda2.o
objs += opkda1.o
objs += opkdmain.o
objs += krome_user_commons.o
objs += krome_all.o
objs += leiden_shielding.o
objs += richings_dissociation_rates.o

exec=test_leiden_shielding

# Include location of KROME
VPATH=../krome-builds/build-HCO-black

default : all

all : test_leiden_shielding test_richings

test_leiden_shielding: $(objs) test_leiden_shielding.o
	$(fc) $(objs) test_leiden_shielding.o -o $(exec) $(switch)

test_richings: exec=test_richings
test_richings: $(objs) test_richings.o
	$(fc) $(objs) test_richings.o -o $(exec) $(switch)

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec)

#rule for f90
%.o:%.f90
	$(fc) $(switch) -c $^ -o $@

#rule for f
%.o:%.f
	$(fc) $(switch) $(nowarn) -c $^ -o $@