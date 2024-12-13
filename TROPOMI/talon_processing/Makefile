# Comment line

EXEC = trop_coloc_thin_exec 
COMP = gfortran
OBJECTS = trop_colocate_thin.o

FFLAGS = -fopenmp -fintrinsic-modules-path /usr/lib64/gfortran/modules -lhdf5_fortran
LDFLAGS = -fopenmp -fintrinsic-modules-path /usr/lib64/gfortran/modules -lhdf5_fortran

$(EXEC): $(OBJECTS)
	$(COMP) $(LDFLAGS) $(OBJECTS) -o $(EXEC)

trop_colocate_thin.o: trop_colocate_thin.f90
	$(COMP) $(FFLAGS) -c trop_colocate_thin.f90

clean:
	rm $(EXEC) $(OBJECTS) *.mod
