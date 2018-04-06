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
#ifdef USE_INTRINSICS
#include "vector_x86.hpp"
#endif
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_hgas.hpp"

class SPH3D : public HelmholtzGas {
public:
    SPH3D() {
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

    void read(FILE * fp) {
        fscanf(fp, "%lld%lld%lf", &this->id, &this->istar, &this->mass);         //  3
        fscanf(fp, "%lf%lf%lf", &this->pos[0], &this->pos[1], &this->pos[2]); //  6
        fscanf(fp, "%lf%lf%lf", &this->vel[0], &this->vel[1], &this->vel[2]); //  9
        fscanf(fp, "%lf%lf%lf", &this->acc[0], &this->acc[1], &this->acc[2]); // 12
        fscanf(fp, "%lf%lf%lf", &this->uene, &this->alph, &this->alphu);      // 15
        fscanf(fp, "%lf%lf%6d", &this->dens, &this->ksr,  &this->np);         // 18
        fscanf(fp, "%lf%lf%lf", &this->vsnd, &this->pres, &this->temp);       // 21
        fscanf(fp, "%lf%lf%lf", &this->divv, &this->rotv, &this->bswt);       // 24
        fscanf(fp, "%lf%lf%lf", &this->pot,  &this->abar, &this->zbar);       // 27
        fscanf(fp, "%lf",       &this->enuc);                                 // 28
        fscanf(fp, "%lf%lf%lf", &this->vsmx, &this->udot, &this->dnuc);       // 31
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) { // 32 -- 44
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        fscanf(fp, "%lf", &this->pot3);
        fscanf(fp, "%lf%lf%lf", &this->tempmax[0], &this->tempmax[1], &this->tempmax[2]);
        fscanf(fp, "%lf", &this->entr);
    }

    void write(FILE * fp = stdout) {
        fprintf(fp, "%8d %2d %+e",  this->id, this->istar, this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->uene, this->alph, this->alphu);
        fprintf(fp, " %+e", this->ksr);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
/////////////// A. Tanikawa adds this 18/03/22
            this->cmps[k] = ((this->cmps[k] != 0.) ? this->cmps[k] : 1e-10);
///////////////
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        fprintf(fp, "\n");
    }
};

namespace HeliumDetonation {
    PS::F64 tmax = 1.0e9;
    PS::F64 tmin = 2.5e8;
    PS::F64 size = 1.0e8;
    PS::F64vec hpos;

    PS::F64 getTemperatureLinear(SPH3D & sph) {
        if(sph.cmps[0] < 0.1) {
            return sph.temp;
        } else {
            PS::F64 temp;
            PS::F64 rsph = sqrt(sph.pos[0] * sph.pos[0] + sph.pos[1] * sph.pos[1]);
            if(rsph < size & sph.pos[2] > 0.) {
                temp = tmax - (tmax - tmin) * rsph / size;
            } else {
                //temp = tmin;
                temp = sph.temp;
            }
            return temp;
        }
    }
/*
    PS::F64 getTemperatureLinear(SPH3D & sph) {
        PS::F64 temp = sph.temp;
        PS::F64vec dx = sph.pos - hpos;
        PS::F64 rsph = sqrt(dx * dx);
        if(rsph < size) {
            temp = tmax - (tmax - tmin) * rsph / size;
            temp = std::max(temp, sph.temp);
        } else {
            temp = tmin;
        }
        return temp;
    }
*/
}

namespace CarbonDetonation {
    PS::F64 tmax = 3.5e9;
    PS::F64 tmin = 1e8;
    PS::F64 size = 0.0;
    PS::F64vec hpos;

    PS::F64 getTemperatureLinear(SPH3D & sph) {
        PS::F64 temp = sph.temp;
        PS::F64vec dx = sph.pos - hpos;
        PS::F64 rsph = sqrt(dx * dx);
        if(rsph < size) {
            temp = tmax - (tmax - tmin) * rsph / size;
            temp = std::max(temp, sph.temp);
        } else {
            temp = tmin;
        }
        temp = (temp > sph.temp) ? temp : sph.temp;
        return temp;
    }
}

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    PS::S64 nptcl;
    char ifile[64];
    char ofile[64];
    PS::S64 flag;
    PS::F64 rmin;
    PS::F64 hemf;

    {
        FILE * fp = fopen(argv[1], "r");
        fscanf(fp, "%lld%s%s", &nptcl, ifile, ofile);
        fscanf(fp, "%lld", &flag);
        printf("nptcl: %8d\n", nptcl);
        printf("ifile: %s\n", ifile);
        printf("ofile: %s\n", ofile);
        if(flag == 0) {
            printf("Helium detonation\n");
            fscanf(fp, "%lf", &HeliumDetonation::size);
            printf("size: %+e\n", HeliumDetonation::size);
            fscanf(fp, "%lf%lf%lf",
                   &HeliumDetonation::hpos[0],
                   &HeliumDetonation::hpos[1],
                   &HeliumDetonation::hpos[2]);
            printf("spot: %+e %+e %+e\n",
                   HeliumDetonation::hpos[0],
                   HeliumDetonation::hpos[1],
                   HeliumDetonation::hpos[2]);
            printf("CAUTION!! fHe < 0.1 particles can not be hotspot!\n");
        } else {
            printf("Carbon detonation\n");
            fscanf(fp, "%lf", &CarbonDetonation::size);
            printf("size: %+e\n", CarbonDetonation::size);
            fscanf(fp, "%lf%lf%lf",
                   &CarbonDetonation::hpos[0],
                   &CarbonDetonation::hpos[1],
                   &CarbonDetonation::hpos[2]);
            printf("spot: %+e %+e %+e\n",
                   CarbonDetonation::hpos[0],
                   CarbonDetonation::hpos[1],
                   CarbonDetonation::hpos[2]);
        }
        fclose(fp);
    }    

    {
        FILE * ifp = fopen(ifile, "r");
        assert(ifp);
        FILE * ofp = fopen(ofile, "w");
        init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);
        if(flag == 0) {
            for(PS::S64 i = 0; i < nptcl; i++) {
                SPH3D sph;
                sph.read(ifp);
                PS::F64 temp = HeliumDetonation::getTemperatureLinear(sph);
                if(temp != sph.temp) {
                    flash_helmholtz_e_(&sph.dens, &temp, sph.cmps.getPointer(), &sph.uene);
                }
                if(sph.uene < 0.) {
                    PS::F64 temp = CodeUnit::MinimumOfTemperature;
                    flash_helmholtz_e_(&sph.dens, &temp, sph.cmps.getPointer(), &sph.uene);
                }
                sph.istar = 0;
                sph.write(ofp);
            }
        } else {
            for(PS::S64 i = 0; i < nptcl; i++) {
                SPH3D sph;
                sph.read(ifp);
                PS::F64 temp = CarbonDetonation::getTemperatureLinear(sph);
                if(temp != sph.temp) {
                    flash_helmholtz_e_(&sph.dens, &temp, sph.cmps.getPointer(), &sph.uene);
                }
                if(sph.uene < 0.) {
                    PS::F64 temp = CodeUnit::MinimumOfTemperature;
                    flash_helmholtz_e_(&sph.dens, &temp, sph.cmps.getPointer(), &sph.uene);
                }
                sph.istar = 0;
                sph.write(ofp);
            }
        }
        fclose(ifp);
        fclose(ofp);
    }

    MPI_Finalize();

    return 0;
}
