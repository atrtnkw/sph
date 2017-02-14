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
    PS::F64 entr;
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

    static PS::F64 getTemperatureImplicitlyWithNse(const SPH sph,
                                                   PS::F64 & dnuc,
                                                   NR::Nucleon & cmps) {
        const PS::S32 nstepmax = 1000;
        const PS::F64 tmin     = CodeUnit::MinimumOfTemperatureNSE;
        const PS::F64 tmax     = CodeUnit::MaximumOfTemperatureNSE;
        PS::S32 nstep   = 0;
        PS::F64 tguess  = sph.temp;
        PS::F64 tresult, tgprev, trtent;
        SPH tsph;
        do {
            tsph       = sph;
            tsph.temp  = tguess;
            tresult    = tguess;
            tsph.dnuc  = CalcNSE::getGeneratedEnergy(tsph.dens, tsph.temp, tsph.cmps.getPointer());
            tsph.uene += tsph.dnuc;
            CalcEquationOfState::getThermodynamicQuantity(tsph.dens, tsph.uene, tsph.cmps,
                                                          tsph.pres, tsph.vsnd, tresult, tsph.entr);
            tgprev = tguess;
            trtent = std::max(tmin, std::min(tmax, tresult));
            PS::F64 rn = getRandomNumber();
            tguess = rn * tguess + (1. - rn) * trtent;
        } while(fabs((tresult - tgprev) / tgprev) > 0.01 && nstep < nstepmax);
        assert(nstep < nstepmax);
        dnuc = tsph.dnuc;
        cmps = tsph.cmps;
        return tresult;
    }

    static PS::F64 getTemperatureImplicitlyWithAprox13(const PS::F64 dt,
                                                       const SPH sph,
                                                       PS::F64 & dnuc,
                                                       NR::Nucleon & cmps) {
        const PS::S32 nstepmax = 1000;
        const PS::F64 tmin     = CodeUnit::MinimumOfTemperatureNSE;
        const PS::F64 tmax     = CodeUnit::MaximumOfTemperatureNSE;
        PS::S32 nstep   = 0;
        PS::F64 tguess  = sph.temp;
        PS::F64 tresult, tgprev, trtent;
        SPH tsph;
        do {
            tsph       = sph;
            tsph.temp  = tguess;
            tresult    = tguess;
            tsph.dnuc  = CalcNRH::getGeneratedEnergy(dt, tsph.dens, tsph.temp, tsph.cmps);
            tsph.uene += tsph.dnuc;
            CalcEquationOfState::getThermodynamicQuantity(tsph.dens, tsph.uene, tsph.cmps,
                                                          tsph.pres, tsph.vsnd, tresult,
                                                          tsph.entr);
            tgprev = tguess;
            trtent = std::max(tmin, std::min(tmax, tresult));
            PS::F64 rn = getRandomNumber();
            tguess = rn * tguess + (1. - rn) * trtent;
        } while(fabs((tresult - tgprev) / tgprev) > 0.01 && nstep < nstepmax);
        assert(nstep < nstepmax);
        dnuc = tsph.dnuc;
        cmps = tsph.cmps;
        return tresult;
    }

    void calcReleasedEnergyConstantTimestep(PS::F64 dt) {
        this->dnuc = CalcNRH::getGeneratedEnergy(dt,
                                                 this->dens,
                                                 this->temp,
                                                 this->cmps);
    }

    void calcReleasedEnergyWithNSE(PS::F64 dt) {
        if(this->temp < CodeUnit::MinimumOfTemperatureNSE) {
            this->dnuc = CalcNRH::getGeneratedEnergy(dt,
                                                     this->dens,
                                                     this->temp,
                                                     this->cmps);
        } else {
            PS::F64 tresult = getTemperatureImplicitlyWithNse(*this, this->dnuc, this->cmps);
        }
    }

/*
    void calcReleasedEnergyHybrid(PS::F64 dt) {
        if(this->temp < CodeUnit::MinimumOfTemperatureNSE
           || this->dens * CodeUnit::UnitOfDensity < 5e7) {
            const PS::F64 dtsm = 1. / (PS::F64)(((PS::S64)1) << 25);
            if(this->temp < CodeUnit::MinimumOfTemperatureNSE || dt <= dtsm) {
                this->dnuc = CalcNRH::getGeneratedEnergy(dt,
                                                         this->dens,
                                                         this->temp,
                                                         this->cmps);
            } else {
                SPH tsph = (*this);
                PS::F64 uinit = tsph.uene;
                for(PS::F64 time = 0.; time < dt; time += dtsm) {
                    tsph.calcReleasedEnergyConstantTimestep(dtsm);
                    tsph.uene += tsph.dnuc;
                    PS::F64 tresult;
                    CalcEquationOfState::getThermodynamicQuantity(tsph.dens, tsph.uene, tsph.cmps,
                                                                  tsph.pres, tsph.vsnd, tresult,
                                                                  tsph.entr);
                }
                this->dnuc = tsph.uene - uinit;
                this->cmps = tsph.cmps;
            }
        } else {
            PS::F64 tresult = getTemperatureImplicitlyWithNse(*this, this->dnuc, this->cmps);
        }
    }
*/

    void calcReleasedEnergyImplicit(PS::F64 dt) {
        if(this->temp < CodeUnit::MinimumOfTemperatureNSE
           || this->dens * CodeUnit::UnitOfDensity < 3e7) {
            this->dnuc = CalcNRH::getGeneratedEnergy(dt,
                                                     this->dens,
                                                     this->temp,
                                                     this->cmps);
        } else {
            PS::F64 tresult = getTemperatureImplicitlyWithAprox13(dt, *this,
                                                                  this->dnuc, this->cmps);
        }
    }

/*
    void calcReleasedEnergySemiImplicit(PS::F64 dt) {
        if(this->temp < CodeUnit::MinimumOfTemperatureNSE) {
            this->dnuc = CalcNRH::getGeneratedEnergy(dt,
                                                     this->dens,
                                                     this->temp,
                                                     this->cmps);
        } else {
            PS::F64     dum_dnuc;
            NR::Nucleon dum_cmps;
            PS::F64 tguess = getTemperatureImplicitlyWithNse(*this, dum_dnuc, dum_cmps);
            this->dnuc = CalcNRH::getGeneratedEnergy(dt, this->dens, tguess, this->cmps);
        }
    }
*/

    void print(FILE *fp,
               PS::F64 time = 0.) {
        fprintf(fp, "%16.16f", time);
        fprintf(fp, " %+e", this->dens * CodeUnit::UnitOfDensity);
        fprintf(fp, " %+e", this->temp);
        fprintf(fp, " %+e", this->uene * CodeUnit::UnitOfEnergy);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) { // 5 -- 17
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        fprintf(fp, "\n");
    }

};

PS::F64 obtainEnergyGivingTemperature(PS::F64 temp,
                                      SPH & sph) {
    PS::F64 dd = sph.dens * CodeUnit::UnitOfDensity;
    /*
    PS::F64 pp = 0.;    
    PS::F64 uu = 0.;
    PS::F64 du = 0.;
    PS::F64 cs = 0.;
    bool eosfail;
    helmeos2_(&temp, &dd, &sph.abar, &sph.zbar, &pp, &uu, &du, &cs, &eosfail);
    */
    PS::F64 uu = 0.;
    flash_helmholtz_e_(&dd, &temp, sph.cmps.getPointer(), &uu);
    PS::F64 uene = uu * CodeUnit::UnitOfEnergyInv;
    return uene;
}

void integrateEnergyByConstantTimestep(FILE *fp,
                                       PS::F64 tend,
                                       PS::F64 dtim,
                                       SPH & sph) {
    fprintf(stderr, "Mode of Constant Timestep\n");
    PS::S64 nstep = 0;
    for(PS::F64 time = 0.; time <= tend; time += dtim) {
        if(nstep % 1024 == 0) {
            sph.print(fp, time);
            fflush(fp);
        }
        sph.calcReleasedEnergyConstantTimestep(dtim);
        sph.uene += sph.dnuc;
        CalcEquationOfState::getThermodynamicQuantity(sph.dens, sph.uene, sph.cmps,
                                                      sph.pres, sph.vsnd, sph.temp,
                                                      sph.entr);
        nstep++;
    }
}

/*
void integrateEnergyWithNse(FILE *fp,
                            PS::F64 tend,
                            PS::F64 dtim,
                            SPH & sph) {
    fprintf(stderr, "Mode with NSE\n");
    for(PS::F64 time = 0.; time <= tend; time += dtim) {
        sph.print(fp, time);
        fflush(fp);
        sph.calcReleasedEnergyWithNSE(dtim);
        sph.uene += sph.dnuc;
        CalcEquationOfState::getThermodynamicQuantity(sph.dens, sph.uene, sph.cmps,
                                                      sph.pres, sph.vsnd, sph.temp,
                                                      sph.entr);
    }
}
*/

/*
void integrateEnergyHybrid(FILE *fp,
                           PS::F64 tend,
                           PS::F64 dtim,
                           SPH & sph) {
    fprintf(stderr, "Mode of Hybrid\n");
    for(PS::F64 time = 0.; time <= tend; time += dtim) {
        sph.print(fp, time);
        fflush(fp);
        sph.calcReleasedEnergyHybrid(dtim);
        sph.uene += sph.dnuc;
        CalcEquationOfState::getThermodynamicQuantity(sph.dens, sph.uene, sph.cmps,
                                                      sph.pres, sph.vsnd, sph.temp,
                                                      sph.entr);
    }
}
*/

void integrateEnergyImplicit(FILE *fp,
                             PS::F64 tend,
                             PS::F64 dtim,
                             PS::S64 dnstep,
                             SPH & sph) {
    fprintf(stderr, "Mode Implicit\n");
    PS::S64 nstep = 0;
    for(PS::F64 time = 0.; time <= tend; time += dtim) {
        if(nstep % dnstep == 0) {
            sph.print(fp, time);
            fflush(fp);
        }
        sph.calcReleasedEnergyImplicit(dtim);
        sph.uene += sph.dnuc;
        CalcEquationOfState::getThermodynamicQuantity(sph.dens, sph.uene, sph.cmps,
                                                      sph.pres, sph.vsnd, sph.temp,
                                                      sph.entr);
        nstep++;
    }
}

/*
void integrateEnergySemiImplicit(FILE *fp,
                                 PS::F64 tend,
                                 PS::F64 dtim,
                                 SPH & sph) {
    fprintf(stderr, "Mode Semi-Implicit\n");
    for(PS::F64 time = 0.; time <= tend; time += dtim) {
        sph.print(fp, time);
        fflush(fp);
        sph.calcReleasedEnergySemiImplicit(dtim);
        sph.uene += sph.dnuc;
        sph.calcAbarZbar();
        CalcEquationOfState::getThermodynamicQuantity(sph.dens, sph.uene, sph.cmps,
                                                      sph.pres, sph.vsnd, sph.temp,
                                                      sph.entr);
    }
}
*/

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    RP::FlagDamping = 0;

    SPH sph;
    PS::S64 te = 0;
    PS::S64 dt = 0;
    PS::S64 dn = 0;
    PS::F64 dd = 0.;
    PS::F64 tt = 0.;
    PS::S64 fg = 0;
    char filename[64];
    FILE * fp = fopen(argv[1], "r");
    //fprintf(stderr, "Input tend=1./(2^xx[s]), dt=1./(2^yy[s])");
    //fprintf(stderr, " temperature[K], density[g/cm^3], FileName, 0/1.\n");
    fscanf(fp, "%lld%lld%lld%lf%lf%s%lld", &te, &dt, &dn, &tt, &dd, filename, &fg);
    if(fg != 0) {
        //fprintf(stderr, "Compositions.\n");
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fscanf(fp, "%lf", &sph.cmps[k]);
        }
    } else {
        sph.cmps[1]  = 0.5;
        sph.cmps[2]  = 0.5;
    }
    fclose(fp);

    init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);

    {
        PS::F64 tend = (1./(PS::F64)(((PS::S64)1) << te)) * CodeUnit::UnitOfTimeInv;
        PS::F64 dtim = (1./(PS::F64)(((PS::S64)1) << dt)) * CodeUnit::UnitOfTimeInv;
        sph.dens     = dd * CodeUnit::UnitOfDensityInv;
        sph.calcAbarZbar();
        sph.uene     = obtainEnergyGivingTemperature(tt, sph);
        CalcEquationOfState::getThermodynamicQuantity(sph.dens, sph.uene, sph.cmps,
                                                      sph.pres, sph.vsnd, sph.temp,
                                                      sph.entr);        
        FILE * fp = fopen(filename, "w");
        PS::F64 st0 = getWallclockTime();
        //integrateEnergyByConstantTimestep(fp, tend, dtim, sph);
        //integrateEnergyWithNse(fp, tend, dtim, sph);
        //integrateEnergyHybrid(fp, tend, dtim, sph);
        integrateEnergyImplicit(fp, tend, dtim, dn, sph);
        //integrateEnergySemiImplicit(fp, tend, dtim, sph);
        PS::F64 st1 = getWallclockTime();
        fclose(fp);
        fprintf(stderr, "WT: %+e\n", st1 - st0);
    }

    MPI_Finalize();

    return 0;
}
