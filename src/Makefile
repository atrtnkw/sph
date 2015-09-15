FDPS_DIR = ../../../fdps-jenever0/src
#FDPS_DIR = ../../fdps/src

PHANTOM_DIR = ${FDPS_DIR}/phantom_grape_x86/G5/newton/libpg5

PS_PATH  = -I${FDPS_DIR}
PS_PATH += -I${PHANTOM_DIR}

CCC = mpicxx
#CCC = CC
#CCC = mpicxx-openmpi-gcc49
#CCC = mpicxx-mpich-gcc49
#CCC = g++

CFLAGS = -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION
#CFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY

CLIBS += -L${PHANTOM_DIR} -lpg5

#########################################
############# AT_FLAGS ##################
#CFLAGS += -DTIMETEST
#
CFLAGS += -DTHERMAL_CONDUCTIVITY
#
CFLAGS += -DENABLE_SIMDX86
#CFLAGS += -DENABLE_SIMDX86_SINGLE
#
CFLAGS += -DWENDLANDC2
#CFLAGS += -DWENDLANDC4
#CFLAGS += -DCUBICSPLINE
#
#CFLAGS += -DUSE_AT1D
CFLAGS += -DUSE_AT2D
#
#CFLAGS += -DGRAVITY
#CFLAGS += -DSYMMETRIZED_GRAVITY
#
#CFLAGS += -DDAMPING
#########################################
#########################################

OBJ = main.o
HDR = hdr_sph.hpp hdr_kernel.hpp hdr_density.hpp hdr_hydro.hpp \
	hdr_gravity.hpp hdr_damp.hpp hdr_time.hpp vector_x86.hpp \
	hdr_vgravity.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

asm: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -S main.cpp

clean:
	rm -f *~
	rm -f test.o*
	rm -f core
	rm -f *.s
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run

change:
	cp hdr_sph.hpp ideal/hdr_sph.hpp
	cp main.cpp ideal/main.cpp
	cp Makefile ideal/Makefile.all
	cp heos/hdr_sph.hpp .
	cp heos/main.cpp .
	cp heos/Makefile.all Makefile
