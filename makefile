# Comment line

COMP = gfortran
OBJECTS = omi_frequency.o \
			count_ai.o \
			synop_time_check.o \
			read_shawn_file.o \
			test_sub.o

omi_exec: $(OBJECTS)
	$(COMP) $(OBJECTS) -o omi_exec 

omi_frequency.o: omi_frequency.f90
	$(COMP) -c omi_frequency.f90

count_ai.o: count_ai.f90
	$(COMP) -c count_ai.f90

synop_time_check.o: synop_time_check.f90
	$(COMP) -c synop_time_check.f90

read_shawn_file.o: read_shawn_file.f90
	$(COMP) -c read_shawn_file.f90

test_sub.o: test_sub.f90
	$(COMP) -c test_sub.f90

clean:
	rm omi_exec
	rm $(OBJECTS)
