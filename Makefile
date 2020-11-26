F90=mpif90
CC=mpicxx
#F90=mpiifort
PETSC_DIR=/home/nvarini/AMGX/petsc-hypre
AMGX_DIR=/home/nvarini/AMGX/v-2.1.x-gnu
AMGXWRAPPER_DIR=/home/nvarini/AMGX/AmgXWrapper
CFLAGS=-I$(AMGX_DIR)/include -I$(AMGXWRAPPER_DIR)/include -I$(PETSC_DIR)/include -g
F90FLAGS=-I$(PETSC_DIR)/include -I  -g -ffree-form -ffree-line-length-none
LDFLAGS=-L$(PETSC_DIR)/lib -lpetsc -L$(AMGX_DIR)/lib -L$(AMGXWRAPPER_DIR)/lib64 -lamgxsh -lAmgXWrapper -lc -lstdc++ -g  
LDFLAGS_INTEL=-L$(PETSC_DIR)/lib -lpetsc -g

          

all:   miniapp_amgx.x miniapp-all.x miniapp-all-mps.x

miniapp-all-mps.x: miniapp-all-mps.cxx petscToCSR.cxx
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
miniapp-all.x: miniapp-all.cxx petscToCSR.cxx
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
miniapp_amgx.x: amgx_mod.o miniapp_amgx.o solver_amgx.o petscToCSR.o
	$(F90) $(F90FLAGS) $(LDFLAGS) amgx_mod.o miniapp_amgx.o solver_amgx.o petscToCSR.o -o $@
miniapp_intel.x:  miniapp_intel.o
	$(F90)  $(LDFLAGS_INTEL) miniapp_intel.o -o $@
amgx_mod.o: amgx_mod.F90
	$(F90) -c $(F90FLAGS) amgx_mod.F90 -o $@
miniapp_amgx.o: miniapp_amgx.F90
	$(F90) -c $(F90FLAGS) miniapp_amgx.F90 -o $@
miniapp_intel.o: miniapp_intel.F90
	$(F90) -c $(F90FLAGS) miniapp_intel.F90 -o $@
solver_amgx.o: solver_amgx.cxx
	$(CC) -c $(CFLAGS) solver_amgx.cxx -o $@
petscToCSR.o: petscToCSR.cxx
	$(CC) -c $(CFLAGS) petscToCSR.cxx -o $@
clean:
	rm *.x *.o
