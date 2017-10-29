###### AT_PLATFORM ########
#use_xc  = yes
use_ofp  = yes
###### AT_EOS #############
use_iso7  = no
use_ideal = no
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
#CFLAGS += -DUSE_INTRINSICS
CFLAGS += -DUSE_HELMHOLTZ
#CFLAGS += -DCOULOMB_CORRECTION
#CFLAGS += -DOLD_FORMAT
#CFLAGS += -DNEW_FORMAT
##############################################
##############################################


ifeq ($(use_xc), yes)
FDPS_DIR = /home/tanikwat/fdps/src
SPH_DIR = /home/tanikwat/git-sph/std
CCC = CC
FCC = ftn
CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CLIBS  += -lstdc++ -lm -L../../flash/ -lflash
endif

ifeq ($(use_ofp), yes)
HOME = /home/x10049
FDPS_DIR = ${HOME}/fdps/src
SPH_DIR = ${HOME}/git-sph/std
CCC = mpiicpc
CNN = mpiicc
FCC = mpiifort
LNK = ${FCC}
CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CFLAGS += -O3 -axMIC-AVX512 #-Wall
CFLAGS += -fimf-precision:low
CLIBS  += -lmpicxx -lstdc++ -nofor_main #-lm 
CLIBS  += -L./flash -lflash
ifeq ($(use_omp), yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -parallel -qopenmp
CFLAGS += -DUSE_XEONPHI
endif
ifeq ($(use_mpi), yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
endif

PS_PATH  = -I${FDPS_DIR}

OBJ  = main.o ${SPH_DIR}/heos/eos.o #${SPH_DIR}/aprox13/aprox13_main.o \
	${SPH_DIR}/aprox13/aprox13_input.o ${SPH_DIR}/aprox13/aprox13.o \
	${SPH_DIR}/nse/nse.o ${SPH_DIR}/nse/nse_main.o
ETC  = Makefile

run: $(OBJ) $(ETC)
	$(FCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

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
