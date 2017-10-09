###### AT_PLATFORM ########
use_xc  = no
use_ofp = yes
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
CFLAGS += -DUSE_HELMHOLTZ
#CFLAGS += -DUSE_INTRINSICS
#CFLAGS += -DFOR_TUBE_TEST
#CFLAGS += -DCOULOMB_CORRECTION
#CFLAGS += -DNON_SYNCHRONIZED_BINARY
#CFLAGS += -DNBODYLIKE
##############################################
##############################################

ifeq ($(use_xc), yes)
HOME_DIR = /home/tanikwat
FDPS_DIR = ${HOME_DIR}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std
CCC = CC
FCC = ftn
CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
ifeq ($(use_omp), yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif
ifeq ($(use_mpi), yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
endif

ifeq ($(use_ofp), yes)
HOME_DIR = /home/x10049
FDPS_DIR = ${HOME_DIR}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std
CCC = mpiicpc
FCC = mpiifort
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

PS_PATH  = -I${FDPS_DIR}

CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CLIBS  += -lstdc++ -lm -L../../flash/ -lflash
CFLAGS += -I${HOME_DIR}/opt/

OBJ  = main.o ${SPH_DIR}/heos/eos.o
ETC  = Makefile

run: $(OBJ) $(ETC)
	$(FCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp $(CLIBS)

clean:
	rm -f *~
	rm -f stdouterr.log
	rm -f core
	rm -f *.s
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run
