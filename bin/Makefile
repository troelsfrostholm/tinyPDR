#default fortran compiler
fc = gfortran

#executable name
exec = tinyPDR

#default libraries
lib = -llapack

#Krome build folder
#KROME_BUILD=/Users/troels/Code/forks/krome/build
#KROME_BUILD=../krome/build_HCO-ice-dust-table
#KROME_BUILD=../krome/build_HCO-noice
KROME_BUILD=../krome/build_HCO-heatcool
#Shielding functions folder
SHIELDING_DIR=../shielding

#Source
SRC=../src

VPATH=$(KROME_BUILD):$(SHIELDING_DIR):$(SRC)

#flags
switchDBG = -fbacktrace -g -O0 -ffpe-trap=zero,overflow,invalid -fbounds-check -ffree-line-length-none -finit-real=nan -cpp

switchDustTable = -DUSE_DUST_TABLE
#switchDustTable=

switchCooling = -DUSE_COOLING
switchHeating = -DUSE_HEATING

#default switch
switch = $(switchDBG) $(switchDustTable) $(switchCooling) $(switchHeating)

#no warning switch (for DLSODES)
nowarnIfort = -nowarn
nowarnGNU = -w

#objects
SRCDIRS := $(subst :, ,$(VPATH))
#files := $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.f90))
EXTENTIONS := .f90 .f
objs := \
$(foreach ext,$(EXTENTIONS),\
	$(foreach dir,$(SRCDIRS),\
		$(patsubst $(dir)/%$(ext),%.o,$(wildcard $(dir)/*$(ext)))))
objs := $(objs)

default: debug

#default target
all: 	$(objs)
	$(fc) $(objs) -o $(exec) $(switch) $(lib)

#ifort full debug target
#debug: switch = $(switchDBG) $(switchSOLVER)
debug: nowarn = $(nowarnGNU)
debug: all

DEPFILE = Makefile.dep
sinclude $(DEPFILE)

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec)

#rule for f90
%.o:%.f90
	$(fc) $(switch) -c $^ -o $@

#rule for f
%.o:%.f
	$(fc) $(switch) $(nowarn) -c $^ -o $@


.PHONY: dependencies
dependencies:
	fortrandep "$(VPATH)" " " > $(DEPFILE)

print-%  : ; @echo $* = $($*)
