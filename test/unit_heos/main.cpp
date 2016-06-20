#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>
#include <mpi.h>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};

#include "particle_simulator.hpp"
#include "hdr_run.hpp"
#include "hdr_time.hpp"
#include "hdr_heos.hpp"
#include "hdr_nuc.hpp"
#include "hdr_nse.hpp"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    RP::FlagDamping = 0;

    init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);
    
    {
        FILE *fp = fopen("result.log", "w");
        PS::F64 tin  = atof(argv[1]);
        NR::Nucleon cmps;
        cmps[0] = 1.0;
//        cmps[1] = cmps[2] = 0.5;
        fprintf(fp, "# din, tin, eout, pout, cout, sout\n");
        for(PS::F64 din = 1e4; din < 1e10; din *= 1.1) {
            PS::F64 eout = 0.0;
            flash_helmholtz_e_(&din, &tin, cmps.getPointer(), &eout);
            PS::F64 ein = eout;
            PS::F64 pout, cout, tout, sout;
            flash_helmholtz_(&din, &ein, &tin, cmps.getPointer(),
                             &pout, &cout, &tout, &sout);
            fprintf(fp, "%+e %+e %+e %+e %+e %+e\n", din, tin, eout, pout, cout, sout);
        }
        fclose(fp);
    }

    MPI_Finalize();

    return 0;
}
