###### AT_PLATFORM ########
use_xc  = no
###########################

ifeq ($(use_xc), yes)
FDPS_DIR = ../../fdps/src
CCC = CC
else
FDPS_DIR = ../../../fdps-jenever0/src
CCC = mpicxx
endif

PHANTOM_DIR = ${FDPS_DIR}/phantom_grape_x86/G5/newton/libpg5

PS_PATH  = -I${FDPS_DIR}
PS_PATH += -I${PHANTOM_DIR}

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
CFLAGS += -DGRAVITY
CFLAGS += -DSYMMETRIZED_GRAVITY
#CFLAGS += -DTHERMAL_CONDUCTIVITY
CFLAGS += -DINTEGRAL_APPROACH_DERIVATIVE
#
CFLAGS += -DENABLE_SIMDX86
#CFLAGS += -DENABLE_SIMDX86_SINGLE
#
CFLAGS += -DWENDLANDC2
#CFLAGS += -DWENDLANDC4
#CFLAGS += -DCUBICSPLINE
#
#CFLAGS += -DUSE_AT1D
#CFLAGS += -DUSE_AT2D
#
#CFLAGS += -DDAMPING
#########################################
#########################################

OBJ = main.o
HDR = hdr_run.hpp vector_x86.hpp hdr_dimension.hpp hdr_kernel.hpp \
	hdr_sph.hpp hdr_igas.hpp hdr_density.hpp hdr_volume.hpp hdr_auxiliary.hpp \
    hdr_hydro.hpp #hdr_gravity.hpp hdr_damp.hpp hdr_time.hpp \
	hdr_vgravity.hpp hdr_gradient.hpp hdr_iad.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

asm: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -S main.cpp

clean:
	rm -f *~
	rm -f *.o[0-9]*
	rm -f core
	rm -f *.s
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run

#change:
#	cp hdr_sph.hpp ideal/hdr_sph.hpp
#	cp main.cpp ideal/main.cpp
#	cp Makefile ideal/Makefile.all
#	cp heos/hdr_sph.hpp .
#	cp heos/main.cpp .
#	cp heos/Makefile.all Makefile
