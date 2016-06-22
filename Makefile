FC = gfortran
FFLAGS = -O3 -g
LAPACK = -llapack

example : example.f90 inverse.o
	$(FC) $(FFLAGS) $(LAPACK) -o example.x example.f90 inverse.o

inverse.o : inverse.f
	$(FC) $(FFLAGS) -c inverse.f
