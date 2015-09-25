###### AT_PLATFORM ########
use_xc  = no
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
CFLAGS += -DENABLE_SIMDX86
#
CFLAGS += -DWENDLANDC2
#CFLAGS += -DWENDLANDC4
#CFLAGS += -DCUBICSPLINE
##############################################
##############################################

ifeq ($(use_xc), yes)
FDPS_DIR = ../../../fdps/src
else
FDPS_DIR = ../../../../fdps-jenever0/src
endif

PS_PATH  = -I${FDPS_DIR}

ifeq ($(use_xc), yes)
CCC = CC
FCC = ftn
else
CCC = g++
FCC = gfortran
ifeq ($(use_mpi), yes)
CCC = mpicxx
FCC = mpif90
endif
endif

CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops
CLIBS  += -lstdc++ -lm

ifeq ($(use_omp), yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif
ifeq ($(use_mpi), yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
ifeq ($(use_xc), yes)
CLIBS  +=
else
CLIBS  += -lmpi_cxx
endif
endif

OBJ = main.o
HDR = hdr_sph.hpp hdr_kernel.hpp hdr_massless.hpp hdr_quantity.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
ifeq ($(use_xc), yes)
	$(CCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)
else
	$(FCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)
endif

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

asm: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -S main.cpp

clean:
	rm -f *~
	rm -f *.s
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
