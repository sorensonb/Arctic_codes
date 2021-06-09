# Comment line

COMP = gfortran
OBJECTS = omi_frequency.o \
			synop_time_check.o \
			test_sub.o

omi_exec: $(OBJECTS)
	$(COMP) $(OBJECTS) -o omi_exec 

omi_frequency.o: omi_frequency.f90
	$(COMP) -c omi_frequency.f90

synop_time_check.o: synop_time_check.f90
	$(COMP) -c synop_time_check.f90

test_sub.o: test_sub.f90
	$(COMP) -c test_sub.f90

clean:
	rm omi_exec
	rm $(OBJECTS)
