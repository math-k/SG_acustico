# Start of the makefile

#Compiler
#COMPILER = f95
COMPILER = ifort

# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
#FCFLAGS = -O3
#FCFLAGS = -O2
#FCFLAGS = -O1 -traceback -debug inline-debug-info
FCFLAGS = -fast -parallel -qopt-report
#FCFLAGS = -fast

#Objects necessary for the main executable
OBJECTS = onda_2d.o module_variables.o
 

onda_2d: $(OBJECTS)
	$(COMPILER) -o onda_2d $(FCFLAGS) $(OBJECTS)

module_variables.mod: module_variables.o module_variables.f90
	$(COMPILER) $(FCFLAGS) -c module_variables.f90

module_variables.o: module_variables.f90
	$(COMPILER) $(FCFLAGS) -c module_variables.f90

onda_2d.o: module_variables.mod onda_2d.f90
	$(COMPILER) $(FCFLAGS) -c onda_2d.f90 

#%.o: %.f90
#	$(COMPILER) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod onda_2d 
	
# End of the makefile
