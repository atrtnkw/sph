
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ######
CFLAGS += -DENABLE_SIMDX86
#CFLAGS += -DWENDLANDC2
#CFLAGS += -DWENDLANDC4
CFLAGS += -DCUBICSPLINE
#CFLAGS += -DUSE_AT1D
#CFLAGS += -DUSE_AT2D
CFLAGS += -DGRAVITY
#CFLAGS += -DWD_DAMPING1
#CFLAGS += -DWD_DAMPING2
#CFLAGS += -DWD_DAMPINGB
######################

FDPS_DIR = ../../../fdps-jenever0/src
#FDPS_DIR = ../../fdps/src

PHANTOM_DIR = ${FDPS_DIR}/phantom_grape_x86/G5/newton/libpg5

PS_PATH  = -I${FDPS_DIR}
PS_PATH += -I${PHANTOM_DIR}

CCC = g++
FCC = gfortran
#CCC = CC
ifeq ($(use_mpi), yes)
CCC = mpicxx
FCC = mpif90
endif

CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
#CFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY
CLIBS  += -L${PHANTOM_DIR} -lpg5
CLIBS  += -lstdc++ -lm
FFLAGS += -I. -I./heos

ifeq ($(use_omp), yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif
ifeq ($(use_mpi), yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CLIBS  += -lmpi_cxx
endif

OBJ = main.o heos/eos.o
HDR = hdr_sph.hpp hdr_eos.hpp hdr_kernel.hpp hdr_density.hpp hdr_hydro.hpp \
	hdr_gravity.hpp hdr_damp.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
#	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)
	$(FCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

heos/eos.o: heos/eos.f
	make -C heos -f Makefile.heos FCC="$(FCC)"

clean:
	rm -f *~
	rm -f test.o*
	rm -f core
	make -C heos -f Makefile.heos clean
distclean:
	make clean
	rm -f *.o
	make -C heos -f Makefile.heos distclean
allclean:
	make distclean
	rm -f run EOS_*.dat
	make -C heos -f Makefile.heos allclean

change:
	cp hdr_sph.hpp heos/hdr_sph.hpp
	cp main.cpp heos/main.cpp
	cp Makefile heos/Makefile.all
	cp ideal/hdr_sph.hpp .
	cp ideal/Makefile.all Makefile
	cp ideal/main.cpp main.cpp
