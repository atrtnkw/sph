#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};

#include "particle_simulator.hpp"
#include "hdr_run.hpp"
#include "hdr_time.hpp"
#include "hdr_heos.hpp"
#include "hdr_nuc.hpp"
#include "hdr_nse.hpp"

inline PS::F64 getRandomNumber() {
    return ((PS::F64) rand() / ((PS::F64)RAND_MAX + 1.));
}

struct SPH {
    PS::F64 dens;
    PS::F64 uene;
    PS::F64 temp;
    PS::F64 pres;
    PS::F64 vsnd;
    PS::F64 abar;
    PS::F64 zbar;
    PS::F64 dnuc;
    NR::Nucleon cmps;
    PS::S64 cnt;

    SPH() {
        dens = 0.;
        uene = 0.;
        temp = 0.;
        pres = 0.;
        vsnd = 0.;
        abar = 0.;
        zbar = 0.;
        dnuc = 0.;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            cmps[k] = 0.;
        }
        cnt  = 0;
    }

    void calcAbarZbar() {
        PS::F64 abarinv = 0.;
        PS::F64 zbar    = 0.;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            abarinv += NR::ainv[k] * this->cmps[k];
            zbar    += NR::zaratio * this->cmps[k];
        }
        this->abar = 1. / abarinv;
        this->zbar = this->abar * zbar;
    }

    void calcReleasedEnergyConstantTimestep(PS::F64 dt) {
        this->dnuc = CalcNRH::getGeneratedEnergy(dt,
                                                 this->dens,
                                                 this->temp,
                                                 this->cmps.getPointer());
    }

    void calcInternalEnergy() {
        PS::F64 dd = this->dens * CodeUnit::UnitOfDensity;
        PS::F64 pp = 0.;    
        PS::F64 uu = 0.;
        PS::F64 du = 0.;
        PS::F64 cs = 0.;
        bool eosfail;
        helmeos2_(&this->temp, &dd, &this->abar, &this->zbar, &pp, &uu, &du, &cs, &eosfail);
        this->uene = uu * CodeUnit::UnitOfEnergyInv;
    }

    void scan(FILE *fp) {
        fscanf(fp, "%lf%lf", &this->temp, &this->dens);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        this->dens *= CodeUnit::UnitOfDensityInv;
    }

    void print(FILE *fp,
               PS::F64 time = 0.) {
        fprintf(fp, "%+e\n", time);
        fprintf(fp, "%+e %+e %+e %+e\n", this->temp,
                this->dens * CodeUnit::UnitOfDensity,
                this->uene * CodeUnit::UnitOfEnergy,
                this->dnuc * CodeUnit::UnitOfEnergy);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) { // 5 -- 17
            fprintf(fp, " %+.8e", this->cmps[k]);
            if(k == 5 || k == 11 || k == NR::NumberOfNucleon - 1) {
                fprintf(fp, "\n");
            }
        }
    }

};

int main(int argc,
         char ** argv) {
    RP::FlagDamping = 0;

    PS::F64 dtime;
    SPH sph;
    FILE *fp = fopen(argv[1], "r");
    fscanf(fp, "%lf", &dtime);
    sph.scan(fp);
    fclose(fp);

    sph.calcAbarZbar();
    sph.calcInternalEnergy();
    sph.calcReleasedEnergyConstantTimestep(dtime);
    
    fp = fopen(argv[2], "w");
    sph.print(fp, dtime);
    fclose(fp);

    return 0;
}
