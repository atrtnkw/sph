###### AT_PLATFORM ########
#use_xc  = yes
use_ofp  = yes
###### AT_EOS #############
use_ideal = no
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
#CFLAGS += -DUSE_INTRINSICS
#CFLAGS += -DFOR_TUBE_TEST
#CFLAGS += -DCOULOMB_CORRECTION
#CFLAGS += -DNON_SYNCHRONIZED_BINARY
#CFLAGS += -DNBODYLIKE
##############################################
##############################################

ifeq ($(use_xc), yes)
HOME_DIR = /home/tanikwat
CCC = CC
CNN = cc
FCC = ftn
LNK = ${CCC}
CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
CLIBS  += -lstdc++ -lm
ifeq ($(use_omp), yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif
ifeq ($(use_mpi), yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
endif

ifeq ($(use_ofp), yes)
HOME_DIR = /home/x10049
CCC = mpiicpc
CNN = mpiicc
FCC = mpiifort
LNK = ${FCC}
CFLAGS += -O3 -axMIC-AVX512 #-Wall
CFLAGS += -fimf-precision:low
CLIBS  += -lmpicxx -lstdc++ -nofor_main #-lm 
ifeq ($(use_omp), yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -parallel -qopenmp
CFLAGS += -DUSE_XEONPHI
endif
ifeq ($(use_mpi), yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
endif

FDPS_DIR = ${HOME_DIR}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std

CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CLIBS  += -lstdc++ -lm -Lflash -lflash

CFLAGS += -DUSE_HELMHOLTZ

PS_PATH  = -I${FDPS_DIR}

OBJ  = main.o ${SPH_DIR}/heos/eos.o
ETC  = Makefile

run: $(OBJ) $(ETC)
ifeq ($(use_xc), yes)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)
endif
ifeq ($(use_ofp), yes)
	$(FCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)
endif

main.o: main.cpp $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp $(CLIBS)

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
