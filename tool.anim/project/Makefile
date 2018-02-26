###### AT_PLATFORM ########
#use_xc  = yes
use_ofp = yes
###### AT_EOS #############
use_iso7  = no
use_ideal = no
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
CFLAGS += -DUSE_HELMHOLTZ
#CFLAGS += -DCOULOMB_CORRECTION
#CFLAGS += -DOLD_FORMAT
#CFLAGS += -DNEW_FORMAT
##############################################
##############################################

ifeq ($(use_xc), yes)
HOME_DIR = /home/tanikwat
FDPS_DIR = ${HOME_DIR}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std
CCC = CC
FCC = ftn
LNK = ${CCC}
CFLAGS += -DUSE_INTRINSICS
ifeq ($(use_omp), yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif
ifeq ($(use_mpi), yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
else
FDPS_DIR = ../../../../fdps-jenever0/src
SPH_DIR = ../../std
CCC = g++
FCC = gfortran
endif

ifeq ($(use_ofp), yes)
HOME_DIR = /home/x10049
FDPS_DIR = ${HOME}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std
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

PS_PATH  = -I${FDPS_DIR}

CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CLIBS  += -lstdc++ -lm -L../../flash/ -lflash
CFLAGS += -I${HOME_DIR}/opt/

OBJ  = main.o ${SPH_DIR}/heos/eos.o #${SPH_DIR}/aprox13/aprox13_main.o \
	${SPH_DIR}/aprox13/aprox13_input.o ${SPH_DIR}/aprox13/aprox13.o \
	${SPH_DIR}/nse/nse.o ${SPH_DIR}/nse/nse_main.o
ETC  = Makefile

run: $(OBJ) $(ETC)
	$(LNK) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

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
