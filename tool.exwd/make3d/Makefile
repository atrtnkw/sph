include ${HOME}/setting/setting.mk
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
##############################################
##############################################

ifeq ($(use_xc), yes)
HOME_DIR = /home/tanikwat
FDPS_DIR = ${HOME_DIR}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std
FLASH_DIR = ${HOME_DIR}/work/git-sph/link/flash
CCC = CC
FCC = ftn
LNK = ${CC}
PS_PATH  = -I${FDPS_DIR}
CFLAGS += -DUSE_INTRINSICS
CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CLIBS  += -lstdc++ -lm -L${FLASH_DIR} -lflash
endif

ifeq ($(use_ofp), yes)
HOME_DIR = ${HOME}
FDPS_DIR = ${HOME_DIR}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std
FLASH_DIR = ${HOME_DIR}/work/git-sph/link/flash
CCC = mpiicpc
FCC = mpiifort
LNK = ${CCC}
PS_PATH  = -I${FDPS_DIR}
CFLAGS += -O3 -axMIC-AVX512 -fimf-precision:low #-Wall
CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CLIBS  += -lmpicxx -lstdc++ -nofor_main -L${FLASH_DIR} -lflash
endif

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
