###### AT_PLATFORM ########
use_xc  = yes
###### AT_EOS #############
use_ideal = no
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
#CFLAGS += -DTIDAL_DISRUPTION_EVENT
##############################################
##############################################

ifeq ($(use_xc), yes)
FDPS_DIR = /home/tanikwat/fdps/src
else
FDPS_DIR = ../../../fdps-jenever0/src
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

CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
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

OBJ  = main.o heos/eos.o flash/libflash.a 	nse/nse_main.o nse/nse.o #\
	aprox13/aprox13_main.o aprox13/aprox13_input.o aprox13/aprox13.o
HDR += hdr_hgas.hpp hdr_heos.hpp hdr_nuc.hpp hdr_util.hpp hdr_nse.hpp
ETC  = Makefile

run: $(OBJ) $(HDR) $(ETC)
ifeq ($(use_xc), yes)
	$(CCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)
else
	$(FCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)
endif

main.o: main.cpp $(HDR) $(LIB) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

heos/eos.o: heos/eos.f
	make -C heos -f Makefile.heos FCC="$(FCC)"

flash/libburner.a:
	make -C flash -f Makefile.local

#aprox13/aprox13_main.o: aprox13/aprox13_main.f90
#	make -C aprox13 -f Makefile.aprox13 aprox13_main.o FCC="$(FCC)"

#aprox13/aprox13_input.o: aprox13/aprox13_input.f90
#	make -C aprox13 -f Makefile.aprox13 aprox13_input.o FCC="$(FCC)"

#aprox13/aprox13.o: aprox13/aprox13.f90
#	make -C aprox13 -f Makefile.aprox13 aprox13.o FCC="$(FCC)"

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
	make -C aprox13 -f Makefile.aprox13 clean
	make -C nse -f Makefile.nse clean
distclean:
	make clean
	rm -f *.o
	make -C heos -f Makefile.heos distclean
	make -C flash -f Makefile.local distclean
	make -C aprox13 -f Makefile.aprox13 distclean	
	make -C nse -f Makefile.nse distclean
allclean:
	make distclean
	rm -f run
	make -C heos -f Makefile.heos allclean
	make -C flash -f Makefile.local allclean
	make -C aprox13 -f Makefile.aprox13 allclean
	make -C nse -f Makefile.nse allclean

endif
