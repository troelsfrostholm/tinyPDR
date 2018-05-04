#default fortran compiler
fc = gfortran

#executable name
exec = tinyPDR

#default libraries
lib = -llapack

#Krome build folder
KROME_BUILD=/Users/troels/Code/forks/krome/build
#KROME_BUILD=../../krome-builds/build-HCO-black

#Shielding functions folder
SHIELDING_DIR=../shielding

#Source
SRC=../src

VPATH=$(KROME_BUILD):$(SHIELDING_DIR):$(SRC)

#flags
switchDBG = -fbacktrace -g -ffpe-trap=zero,overflow,invalid -fbounds-check -ffree-line-length-none
switchSOLVER = -D$(SOLVER) -cpp

#default switch
switch = $(switchDBG) $(switchSOLVER)

#no warning switch (for DLSODES)
nowarnIfort = -nowarn
nowarnGNU = -w

#objects
objs = opkda2.o
objs += opkda1.o
objs += opkdmain.o
objs += krome_user_commons.o
objs += krome_all.o
objs += leiden_shielding.o
objs += richings_dissociation_rates.o
objs += main.o

default: debug

#default target
all: 	$(objs)
	$(fc) $(objs) -o $(exec) $(switch) $(lib)

#ifort full debug target
#debug: switch = $(switchDBG) $(switchSOLVER)
debug: nowarn = $(nowarnGNU)
debug: all

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec)

#rule for f90
%.o:%.f90
	$(fc) $(switch) -c $^ -o $@

#rule for f
%.o:%.f
	$(fc) $(switch) $(nowarn) -c $^ -o $@
