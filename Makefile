FC = gfortran
FFLAGS = -O3 -g
LAPACK = -llapack

SRC = inverse.f90 example.f90
OBJ = $(SRC:%.f90=%.o)

example.x: $(OBJ)
	$(FC) $(FFLAGS) $^ -o $@ $(LAPACK)

$(OBJ): %.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm *.o *.x *.mod