#FDPS_DIR = /nfshome/tanikawa/work3/fdps-arma/src
FDPS_DIR = ../../../fdps-jenever0/src
PHANTOM_DIR = ${FDPS_DIR}/phantom_grape_x86/G5/newton/libpg5

PS_PATH  = -I${FDPS_DIR}
PS_PATH += -I${PHANTOM_DIR}

CCC = mpicxx
#CCC = CC
#CCC = mpicxx-openmpi-gcc49
#CCC = mpicxx-mpich-gcc49
#CCC = g++

CFLAGS = -O3 -ffast-math -funroll-loops #-Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION
CFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY

CLIBS += -L${PHANTOM_DIR} -lpg5

###### AT_FLAGS ######
#CFLAGS += -DWENDLANDC2
CFLAGS += -DCUBICSPLINE
#CFLAGS += -DUSE_AT1D
#CFLAGS += -DUSE_AT2D
#CFLAGS += -DGRAVITY
#CFLAGS += -DDAMPING
######################

OBJ = main.o
HDR = hdr_sph.hpp hdr_kernel.hpp hdr_density.hpp hdr_hydro.hpp \
	hdr_gravity.hpp hdr_damp.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

clean:
	rm -f *~
#	rm -f snap/*
	rm -f test.o*
	rm -f core
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run
