FC= gfortran
FFLAGS = -O3 --free-form -fdefault-real-8 -fimplicit-none

SRC=src/

MODSRC = $(addprefix $(SRC), threedderivatives.f90 onedpoissonsolver.f90 twodpoissonsolver.f90 threedpoissonsolver.f90)
MODULES = $(MODSRC:.f90=.mod)
TEST1SRC = $(addprefix $(SRC),test1dpoisson.f90)
TEST2SRC = $(addprefix $(SRC),test2dpoisson.f90)
TEST3SRC = $(addprefix $(SRC),test3dpoisson.f90)
TESTDERSRC = $(addprefix $(SRC),testthreedder.f90)
OBJECTS = $(TEST1SRC:.f90=.o) $(TEST2SRC:.f90=.o) $(TEST3SRC:.f90=.o) $(TESTDERSRC:.f90=.o)

default : all

all : test1 test2 test3 testder

clean :
	rm $(OBJECTS)
	rm $(MODULES)

$(MODULES) : %.mod:%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJECTS) : %.o:%.f90 $(MODULES)
	$(FC) $(FFLAGS) -c $< -o $@

$(PROGRAM) : $(OBJECTS) $(MODULES)
	$(FC) $(FFLAGS) $(MODULES) $(OBJECTS) -o $@
	mv $@ .

test1 : $(OBJECTS) $(MODULES)
	$(FC) $(FFLAGS) $(MODULES) $(TEST1SRC:.f90=.o) -o test1dpoison

test2 : $(OBJECTS) $(MODULES)
	$(FC) $(FFLAGS) $(MODULES) $(TEST2SRC:.f90=.o) -o test2dpoison

test3 : $(OBJECTS) $(MODULES)
	$(FC) $(FFLAGS) $(MODULES) $(TEST3SRC:.f90=.o) -o test3dpoison

testder : $(OBJECTS) $(MODULES)
	$(FC) $(FFLAGS) $(MODULES) $(TESTDERSRC:.f90=.o) -o test3dder
