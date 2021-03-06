###### AT_PLATFORM ########
use_xc  = no
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
CFLAGS += -DTHERMAL_CONDUCTIVITY
CFLAGS += -DSYMMETRIZED_GRAVITY
#CFLAGS += -DNUCLEAR_REACTION
#
CFLAGS += -DENABLE_SIMDX86
#CFLAGS += -DENABLE_SIMDX86_SINGLE
#
CFLAGS += -DWENDLANDC2
#CFLAGS += -DWENDLANDC4
#CFLAGS += -DCUBICSPLINE
#
#CFLAGS += -DWD_DAMPING1
#CFLAGS += -DWD_DAMPING2
#CFLAGS += -DWD_DAMPINGB
##############################################
##############################################

ifeq ($(use_xc), yes)
FDPS_DIR = ../../fdps/src
else
FDPS_DIR = ../../../fdps-jenever0/src
endif

PHANTOM_DIR = ${FDPS_DIR}/phantom_grape_x86/G5/newton/libpg5

PS_PATH  = -I${FDPS_DIR}
PS_PATH += -I${PHANTOM_DIR}

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
ifeq ($(use_xc), yes)
CLIBS  +=
else
CLIBS  += -lmpi_cxx
endif
endif

OBJ = main.o heos/eos.o aprox13/aprox13_main.o aprox13/aprox13_input.o \
	aprox13/aprox13.o
HDR = hdr_sph.hpp hdr_eos.hpp hdr_kernel.hpp hdr_density.hpp hdr_hydro.hpp \
	hdr_gravity.hpp hdr_damp.hpp hdr_time.hpp hdr_massless.hpp hdr_nuc.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
ifeq ($(use_xc), yes)
	$(CCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)
else
	$(FCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)
endif

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

heos/eos.o: heos/eos.f
	make -C heos -f Makefile.heos FCC="$(FCC)"

aprox13/aprox13_main.o: aprox13/aprox13_main.f90
	make -C aprox13 -f Makefile.aprox13 aprox13_main.o FCC="$(FCC)"

aprox13/aprox13_input.o: aprox13/aprox13_input.f90
	make -C aprox13 -f Makefile.aprox13 aprox13_input.o FCC="$(FCC)"

aprox13/aprox13.o: aprox13/aprox13.f90
	make -C aprox13 -f Makefile.aprox13 aprox13.o FCC="$(FCC)"

asm: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -S main.cpp

clean:
	rm -f *~
	rm -f *.o[0-9]*
	rm -f core
	rm -f *.s
	make -C heos -f Makefile.heos clean
	make -C aprox13 -f Makefile.aprox13 clean
distclean:
	make clean
	rm -f *.o
	make -C aprox13 -f Makefile.aprox13 distclean
allclean:
	make distclean
	rm -f run EOS_*.dat
	make -C aprox13 -f Makefile.aprox13 allclean

change:
	cp hdr_sph.hpp heos/hdr_sph.hpp
	cp main.cpp heos/main.cpp
	cp Makefile heos/Makefile.all
	cp ideal/hdr_sph.hpp .
	cp ideal/Makefile.all Makefile
	cp ideal/main.cpp main.cpp
