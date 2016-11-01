#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>
#include <mpi.h>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};

#include "particle_simulator.hpp"
#include "hdr_time.hpp"
#include "hdr_run.hpp"
#include "vector_x86.hpp"
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_hgas.hpp"

namespace HotSpot {
    PS::F64 tmax;
    PS::F64 tmin;
    PS::F64 width;
}

class SPH1D : public HelmholtzGas {
public:
    SPH1D() {
        this->id    = 0;
        this->istar = 0;
        this->mass  = 0.;
        this->pos   = 0.;
        this->vel   = 0.;
        this->uene  = 0.;
        this->alph  = 0.;
        this->alphu = 0.;
        this->ksr   = 0.;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            this->cmps[k] = 0.;
        }        
    }

    void print(FILE * fp = stdout) {
        fprintf(fp, "%8d %2d %+e",  this->id, this->istar, this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->uene, this->alph, this->alphu);
        fprintf(fp, " %+e", this->ksr);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        fprintf(fp, "\n");
    }
};

PS::F64 getTemperatureTopHat(PS::F64 x) {
    using namespace HotSpot;
    PS::F64 xmax = width * 0.5;

    PS::F64 xabs = fabs(x);
    PS::F64 temp = 0.;
    if(xabs < xmax) {
        temp = tmax;
    } else {
        temp = tmin;
    }
    return temp;
}

PS::F64 getTemperatureLinear(PS::F64 x) {
    /*
    const PS::F64 tmax  = 3e9;
    const PS::F64 tmin  = 1e7;
    const PS::F64 width = 1e8;
    const PS::F64 xmax  = width * 0.5;
    */
    using namespace HotSpot;
    PS::F64 xmax = width * 0.5;

    PS::F64 xabs = fabs(x);
    PS::F64 temp = 0.;
    if(xabs < xmax) {
        temp = tmax - (tmax - tmin) * xabs / xmax;
    } else {
        temp = tmin;
    }
    return temp;
}

int main(int argc, char ** argv) {
//    PS::Initialize(argc, argv);
    MPI_Init(&argc, &argv);

    init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);

    PS::F64 dn_a;
    PS::F64 dn_b;
    PS::F64 vz_c;
    PS::S64 nhlf;
    PS::F64 zmax;
    PS::F64 temp;
    PS::F64 leng;
    NR::Nucleon cmps;
    PS::F64 fhe, fcb, fox, fne, fmg;
    char filetype[1024];

    {
        FILE *fin = fopen(argv[1], "r");
        if(fin == NULL) {
            fprintf(stderr, "Not found file %s\n", argv[1]);
            MPI_Finalize();
            exit(0);
        }
        fscanf(fin, "%lf%lf%lf", &dn_a, &dn_b, &vz_c);
        fscanf(fin, "%lf%lf%lf%lf%lf", &fhe, &fcb, &fox, &fne, &fmg);
        fscanf(fin, "%lld%lf%lf%lf", &nhlf, &temp, &zmax, &leng);
        fscanf(fin, "%s", filetype);
        fclose(fin);
        cmps[0] = fhe;
        cmps[1] = fcb;
        cmps[2] = fox;
        cmps[3] = fne;
        cmps[4] = fmg;
    }

    PS::F64 *pz = new PS::F64[nhlf*2];
    PS::F64 *vz = new PS::F64[nhlf*2];
    PS::F64 *ue = new PS::F64[nhlf*2];
    PS::F64 *hs = new PS::F64[nhlf*2];

    PS::F64 mtot = -dn_a / 3. * zmax * zmax * zmax + dn_b * zmax;
    PS::F64 m    = mtot / nhlf;
    for(PS::S64 i = 0; i < nhlf; i++) {
        PS::F64 zp   = 0.;
        PS::F64 dens = - dn_a * zp * zp + dn_b;
        if(i != 0) {
            zp    = pz[i-1];
            dens  = - dn_a * zp * zp + dn_b;
            pz[i] = zp + m / dens;
        }
        vz[i] = - vz_c * pz[i];
        flash_helmholtz_e_(&dens, &temp, cmps.getPointer(), &ue[i]);
        hs[i] = 5. * m / dens;
    }

    char filename[1024];
    sprintf(filename, "%s.data", filetype);
    FILE * fp = fopen(filename, "w");
    for(PS::S64 i = nhlf - 1; i >= 0; i--) {
        if(i == 0) {
            continue;
        }
        fprintf(fp, " %8d %2d %+e", nhlf - 1 - i, 0, m);
        fprintf(fp, " %+e %+e %+e", -pz[i], 0.0, 0.0);
        fprintf(fp, " %+e %+e %+e", -vz[i], 0.0, 0.0);
        fprintf(fp, " %+e %+e %+e",  ue[i], 1.0, 0.0);
        fprintf(fp, " %+e", hs[i]);
        /*
        fprintf(fp, " 0.0 0.6 0.35 0.05");
        for(PS::S64 k = 0; k < 9; k++) {
            fprintf(fp, " 0.0");
        }
        */
        fprintf(fp, " %+.3e %+.3e %+.3e %+.3e %+.3e", cmps[0], cmps[1], cmps[2], cmps[3], cmps[4]);
        for(PS::S64 k = 0; k < 8; k++) {
            fprintf(fp, " 0.0");
        }
        fprintf(fp, "\n");
    }
    for(PS::S64 i = 0; i < nhlf; i++) {
        fprintf(fp, " %8d %2d %+e", i + nhlf - 1, 0, m);
        fprintf(fp, " %+e %+e %+e", pz[i], 0.0, 0.0);
        fprintf(fp, " %+e %+e %+e", vz[i], 0.0, 0.0);
        fprintf(fp, " %+e %+e %+e", ue[i], 1.0, 0.0);
        fprintf(fp, " %+e", hs[i]);
        /*
        fprintf(fp, " 0.0 0.6 0.35 0.05");
        for(PS::S64 k = 0; k < 9; k++) {
            fprintf(fp, " 0.0");
        }
        */
        fprintf(fp, " %+.3e %+.3e %+.3e %+.3e %+.3e", cmps[0], cmps[1], cmps[2], cmps[3], cmps[4]);
        for(PS::S64 k = 0; k < 8; k++) {
            fprintf(fp, " 0.0");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    sprintf(filename, "%s.oned", filetype);
    fp = fopen(filename, "w");
    fprintf(fp, "%e\n", leng);
    fclose(fp);

    MPI_Finalize();

    return 0;
}
