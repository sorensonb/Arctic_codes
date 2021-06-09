# Comment line

COMP = gfortran
OBJECTS = omi_frequency.o \
			mie_calc.o

omi_exec: $(OBJECTS)
	$(COMP) $(OBJECTS) -o omi_exec 

mie_calc.o: mie_calc.f90
	$(COMP) -c mie_calc.f90

omi_frequency.o: omi_frequency.f90
	$(COMP) -c omi_frequency.f90

clean:
	rm omi_exec
	rm $(OBJECTS)
