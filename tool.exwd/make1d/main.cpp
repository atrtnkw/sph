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

    PS::S64 tnptcl;
    PS::S64 typeperb;
    PS::F64 tlength;
    PS::F64 tdensity;
    PS::F64 talph;
    PS::F64 the4, tc12, to16;
    char filetype[64];

    {
        FILE *fin = fopen(argv[1], "r");
        fscanf(fin, "%lf%lf%lf",    &tlength, &tdensity, &talph);
        fscanf(fin, "%lf%lf%lf",    &the4, &tc12, &to16);
        fscanf(fin, "%lf%lf%lf", &HotSpot::tmax, &HotSpot::tmin, &HotSpot::width);
        fscanf(fin, "%lld", &typeperb);
        fscanf(fin, "%lld", &tnptcl);
        fscanf(fin, "%s", filetype);
        fclose(fin);
    }    

    const PS::S64 nptcl   = tnptcl;
    const PS::F64 length  = tlength;
    const PS::F64 density = tdensity;
    const PS::F64 alph    = talph;

    PS::F64 (*func)(PS::F64);
    if(typeperb == 0) {
        func = getTemperatureLinear;
    } else {
        func = getTemperatureTopHat;
    }

    NR::Nucleon cmps;
    cmps[0] = the4;
    cmps[1] = tc12;
    cmps[2] = to16;

    init_flash_helmholtz_();

    char filename[64];
    sprintf(filename, "%s.data", filetype);
    FILE * fp = fopen(filename, "w");
    for(PS::S64 i = 0; i < nptcl; i++) {
        SPH1D sph;
        sph.id     = i;
        sph.mass   = length * density / nptcl;
        sph.dens   = density;
        sph.pos[0] = i  * length / nptcl - 0.5 * length;
        sph.alph   = alph;
        sph.ksr    = 5. * length / nptcl;
        sph.cmps   = cmps;
        {
            PS::F64 pp, du, cs;
            bool eosfail;
            PS::F64 temp = func(sph.pos[0]);
            flash_helmholtz_e_(&sph.dens, &temp, sph.cmps.getPointer(), &sph.uene);
        }
        sph.print(fp);
    }
    fclose(fp);

    sprintf(filename, "%s.oned", filetype);
    fp = fopen(filename, "w");
    fprintf(fp, "%e\n", length);
    fclose(fp);

//    PS::Finalize();
    MPI_Finalize();

    return 0;
}
