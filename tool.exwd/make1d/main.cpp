#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>

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
    PS::S64 tnptcl;
    PS::F64 tlength;
    PS::F64 tdensity;
    char filetype[64];

    {
        FILE *fin = fopen(argv[1], "r");
        fscanf(fin, "%lf%lf",    &tlength, &tdensity);
        fscanf(fin, "%lf%lf%lf", &HotSpot::tmax, &HotSpot::tmin, &HotSpot::width);
        fscanf(fin, "%lld", &tnptcl);
        fscanf(fin, "%s", filetype);
        fclose(fin);
    }    

    const PS::S64 nptcl   = tnptcl;
    const PS::F64 length  = tlength;
    const PS::F64 density = tdensity;

    PS::F64 (*func)(PS::F64) = getTemperatureLinear;

    NR::Nucleon cmps;
    cmps[1] = cmps[2] = 0.5;

    char filename[64];
    sprintf(filename, "%s.data", filetype);
    FILE * fp = fopen(filename, "w");
    for(PS::S64 i = 0; i < nptcl; i++) {
        SPH1D sph;
        sph.id     = i;
        sph.mass   = length * density / nptcl;
        sph.dens   = density;
        sph.pos[0] = i  * length / nptcl - 0.5 * length;
        sph.alph   = 1.;
        sph.ksr    = 5. * length / nptcl;
        sph.cmps   = cmps;
        {
            PS::F64 pp, du, cs;
            bool eosfail;
            PS::F64 temp = func(sph.pos[0]);
            sph.calcAbarZbar();
            helmeos2_(&temp, &sph.dens, &sph.abar, &sph.zbar, &pp, &sph.uene, &du, &cs, &eosfail);
        }
        sph.print(fp);
    }
    fclose(fp);

    sprintf(filename, "%s.oned", filetype);
    fp = fopen(filename, "w");
    fprintf(fp, "%e\n", length);
    fclose(fp);

    return 0;
}
