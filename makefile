# Comment line

COMP = gfortran
OBJECTS = omi_frequency.o

omi_exec: $(OBJECTS)
	$(COMP) $(OBJECTS) -o omi_exec 

omi_frequency.o: omi_frequency.f90
	$(COMP) -c omi_frequency.f90

clean:
	rm omi_exec
	rm $(OBJECTS)
