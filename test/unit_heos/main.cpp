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

    PS::S64 wdtype = atoi(argv[1]);

    init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);
    
    {
        FILE *fp = fopen("result.log", "w");
        NR::Nucleon cmps;
        PS::S64 ntin;
        PS::F64 *atin;
        if(wdtype == 0) {
            fprintf(fp, "# WD type: HeWD\n");
            cmps[0] = 1.0;
            ntin = 4;
            atin = (PS::F64 *)malloc(ntin*sizeof(PS::F64));
            atin[0] = 1e6;
            atin[1] = 1e7;
            atin[2] = 1e8;
            atin[3] = 3e8;
        } else if (wdtype == 1) {
            fprintf(fp, "# WD type: COWD\n");
            cmps[1] = cmps[2] = 0.5;
            ntin = 5;
            atin = (PS::F64 *)malloc(ntin*sizeof(PS::F64));
            atin[0] = 1e6;
            atin[1] = 1e7;
            atin[2] = 1e8;
            atin[3] = 1e9;
            atin[4] = 3e9;
        } else {
            fprintf(fp, "# WD type: ONeMgWD\n");
            cmps[2] = 0.6;
            cmps[3] = 0.35;
            cmps[4] = 0.05;
            ntin = 5;
            atin = (PS::F64 *)malloc(ntin*sizeof(PS::F64));
            atin[0] = 1e6;
            atin[1] = 1e7;
            atin[2] = 1e8;
            atin[3] = 1e9;
            atin[4] = 3e9;
        }
        fprintf(fp, "# din, tin, eout, pout, cout, sout\n");
        for(PS::S64 itin = 0; itin < ntin; itin++) {
            PS::F64 tin = atin[itin];
            for(PS::F64 din = 1e4; din < 1e10; din *= 1.1) {
                PS::F64 eout = 0.0;
                flash_helmholtz_e_(&din, &tin, cmps.getPointer(), &eout);
                PS::F64 ein = eout;
                PS::F64 pout, cout, tout, sout;
                flash_helmholtz_(&din, &ein, &tin, cmps.getPointer(),
                                 &pout, &cout, &tout, &sout);
                fprintf(fp, "%+e %+e %+e %+e %+e %+e\n", din, tin, eout, pout, cout, sout);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    MPI_Finalize();

    return 0;
}
