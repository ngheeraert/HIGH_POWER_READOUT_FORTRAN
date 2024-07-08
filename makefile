F90 = ifort
#F90 = gfortran
LLAPACK  = -mkl
#LLAPACK  = -llapack
DEFS = -DDP 

OPTFLAGS = -O2 #-fall-intrinsics# -Wl#,-rpath,${MKLROOT}/lib #-fall-intrinsics -O2 -fcheck=all #-Wall -Wtabs
#OPTFLAGS = -O3 #-parallel#-fall-intrinsics -O2 -fcheck=all #-Wall -Wtabs
#OPTFLAGS = #-fall-intrinsics -O2 -fcheck=all #-Wall -Wtabs
FFLAGS  = -cpp $(DEFS) $(OPTFLAGS) $(MPFUN) $(MKLINCLUDE)
LDFLAGS = $(OPTFLAGS) $(LLAPACK)

###### Section for developers

.SUFFIXES: 
.SUFFIXES: .o .f90 .mod

## main modules and objects ##
EXE  = mpol
OBJC = typedefs.o lapackblas.o systm.o dynamics.o main.o 
MODF = typedefs.mod lapackblas.mod systm.mod dynamics.mod 


.f90.mod:
	$(F90) -c $(FFLAGS) $<

.f90.o:
	$(F90) -c $(FFLAGS) $<

mpol: $(MODC) $(OBJC)
	$(F90) -o $@ $(OBJC) $(LDFLAGS)

clean:
	rm -f *.o *.mod *~
 
