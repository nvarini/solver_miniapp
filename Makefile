F90=mpif90
CC=mpicxx
PETSC_DIR=/home/nvarini/miniapp_amgx_new/petsc-amgx-install-gnu
AMGX_DIR=/home/nvarini/miniapp_amgx_new/amgx-install
AMGXWRAPPER_DIR=/home/nvarini/miniapp_amgx_new/amgxwrapper-install
CFLAGS=-I$(AMGX_DIR)/include -I$(AMGXWRAPPER_DIR)/include -I$(PETSC_DIR)/include -g
F90FLAGS=-I$(PETSC_DIR)/include  -g -ffree-form -ffree-line-length-none
LDFLAGS=-L$(PETSC_DIR)/lib -lpetsc -L$(AMGX_DIR)/lib -L$(AMGXWRAPPER_DIR)/lib64 -lamgxsh -lAmgXWrapper -lc -lstdc++

          

all:   miniapp_amgx.x
miniapp_amgx.x: miniapp_amgx.o solver_amgx.o petscToCSR.o
	$(F90) $(F90FLAGS) $(LDFLAGS) miniapp_amgx.o solver_amgx.o petscToCSR.o -o $@
miniapp_amgx.o: miniapp_amgx.F90
	$(F90) -c $(F90FLAGS) miniapp_amgx.F90 -o $@
solver_amgx.o: solver_amgx.cxx
	$(CC) -c $(CFLAGS) solver_amgx.cxx -o $@
petscToCSR.o: petscToCSR.cxx
	$(CC) -c $(CFLAGS) petscToCSR.cxx -o $@
clean:
	rm *.x *.o
