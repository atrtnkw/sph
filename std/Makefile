###### AT_PLATFORM ########
use_xc  = yes
#use_ofp  = yes
###### AT_EOS #############
use_ideal = no
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
CFLAGS += -DUSE_INTRINSICS
#CFLAGS += -DFOR_TUBE_TEST
#CFLAGS += -DCOULOMB_CORRECTION
#CFLAGS += -DNON_SYNCHRONIZED_BINARY
#CFLAGS += -DNBODYLIKE
##############################################
##############################################

ifeq ($(use_xc), yes)
FDPS_DIR = /home/tanikwat/fdps/src
CCC = CC
CNN = cc
FCC = ftn
LNK = ${CCC}
#CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
CFLAGS += -O3 -march=skylake-avx512 -ffast-math -funroll-loops
CLIBS  += -lstdc++ -lm
ifeq ($(use_omp), yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif
ifeq ($(use_mpi), yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
endif

ifeq ($(use_ofp), yes)
FDPS_DIR = /home/x10049/fdps/src
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

HDR = hdr_run.hpp hdr_time.hpp vector_x86.hpp hdr_dimension.hpp hdr_kernel.hpp \
	hdr_sph.hpp hdr_msls.hpp hdr_density.hpp hdr_hydro.hpp hdr_gravity.hpp \
	hdr_bhns.hpp

ifeq ($(use_ideal), yes)

CFLAGS += -DUSE_IDEAL
OBJ  = main.o
HDR += hdr_igas.hpp
ETC  = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

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

else

FFLAGS += -I. -I./heos
CFLAGS += -DUSE_HELMHOLTZ
CLIBS  += -L./flash -lflash

OBJ  = main.o heos/eos.o flash/libflash.a 	nse/nse_main.o nse/nse.o
HDR += hdr_hgas.hpp hdr_heos.hpp hdr_nuc.hpp hdr_util.hpp hdr_nse.hpp
ETC  = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(LNK) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(HDR) $(LIB) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

heos/eos.o: heos/eos.f
	make -C heos -f Makefile.heos FCC="$(FCC)"

flash/libflash.a:
	make -C flash -f Makefile.local FCC="$(FCC)" CNN="$(CNN)"

nse/nse_main.o: nse/nse_main.f
	make -C nse -f Makefile.nse nse_main.o FCC="$(FCC)"

nse/nse.o: nse/nse.f
	make -C nse -f Makefile.nse nse.o FCC="$(FCC)"

asm: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -S main.cpp

clean:
	rm -f *~
	rm -f *.o[0-9]*
	rm -f core
	rm -f *.s
	make -C heos -f Makefile.heos clean
	make -C flash -f Makefile.local clean
	make -C nse -f Makefile.nse clean
distclean:
	make clean
	rm -f *.o
	make -C heos -f Makefile.heos distclean
	make -C flash -f Makefile.local distclean
	make -C nse -f Makefile.nse distclean
allclean:
	make distclean
	rm -f run
	make -C heos -f Makefile.heos allclean
	make -C flash -f Makefile.local allclean
	make -C nse -f Makefile.nse allclean

endif
