#pragma once

#include "hdr_heos.hpp"
#include "hdr_nuc.hpp"
#include "hdr_util.hpp"

static PS::U64 convertF64ToU64(PS::F64 val) {
    union converter {
        PS::F64 f;
        PS::U64 u;
    };
    union converter var;
        var.f = val;
        return var.u;
}
static PS::F64 convertU64ToF64(PS::U64 val) {
    union converter {
        PS::F64 f;
        PS::U64 u;
    };
    union converter var;
    var.u = val;
    return var.f;
}

class HelmholtzGas : public SPH {
public:
    PS::S64 istar;
    PS::F64 temp;
    PS::F64 cmps[NR::NumberOfNucleon];
    PS::S64 cnteos;    
    PS::F64 dnuc;
    PS::F64 enuc;

    void readAscii(FILE * fp) {
        using namespace CodeUnit;
        fscanf(fp, "%lld%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->istar, &this->mass,        // 3
               &this->pos[0], &this->pos[1], &this->pos[2], // 6
               &this->vel[0], &this->vel[1], &this->vel[2], // 9
               &this->uene,   &this->alph,   &this->alphu,  // 12
               &this->ksr);                                 // 13
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {  // 26
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        this->mass *= UnitOfMassInv;
        this->pos  *= UnitOfLengthInv;
        this->vel  *= UnitOfVelocityInv;
        this->uene *= UnitOfEnergyInv;
        this->ksr  *= UnitOfLengthInv;
    }
    
    void writeAscii(FILE *fp) const {
        using namespace CodeUnit;
        PS::F64    tmass = this->mass * UnitOfMass;
        PS::F64vec tpos  = this->pos  * UnitOfLength;
        PS::F64vec tvel  = this->vel  * UnitOfVelocity;
        PS::F64vec tacc  = this->acc  * UnitOfAcceleration;
        PS::F64    tuene = this->uene * UnitOfEnergy;
        PS::F64    tksr  = this->ksr  * UnitOfLength;
        PS::F64    tdens = this->dens * UnitOfDensity;
        PS::F64    tvsnd = this->vsnd * UnitOfVelocity;
        PS::F64    tpres = this->pres * UnitOfPressure;
        PS::F64    tdivv = this->divv * UnitOfTimeInv; // divv [s^-1]
        PS::F64    trotv = this->rotv * UnitOfTimeInv; // rotv [s^-1]
        PS::F64    tpot  = this->pot  * UnitOfEnergy;
        PS::F64    tenuc = this->enuc * UnitOfEnergy;
        PS::F64    tvsmx = this->vsmx * UnitOfVelocity;
        PS::F64    tudot = this->udot * UnitOfEnergy * UnitOfTimeInv;
        PS::F64    tdnuc = this->dnuc * UnitOfEnergy;
        fprintf(fp, "%6d %2d %+e", this->id, this->istar, tmass);    //  3
        fprintf(fp, " %+e %+e %+e", tpos[0], tpos[1], tpos[2]);      //  6
        fprintf(fp, " %+e %+e %+e", tvel[0], tvel[1], tvel[2]);      //  9
        fprintf(fp, " %+e %+e %+e", tacc[0], tacc[1], tacc[2]);      // 12
        fprintf(fp, " %+e %+e %+e", tuene, this->alph, this->alphu); // 15
        fprintf(fp, " %+e %+e %6d", tdens, tksr,  this->np);         // 18
        fprintf(fp, " %+e %+e %+e", tvsnd, tpres, this->temp);       // 21
        fprintf(fp, " %+e %+e %+e", tdivv, trotv, this->bswt);       // 24
        fprintf(fp, " %+e %+e %+e", tpot, this->abar, this->zbar);   // 27
        fprintf(fp, " %+e",         tenuc);                          // 28
        fprintf(fp, " %+e %+e %+e", tvsmx, tudot, tdnuc);            // 31
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) { // 32 -- 44
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        if(RP::FlagDamping == 2) {
            PS::F64vec tomg = RP::RotationalVelocity * UnitOfTimeInv;
            PS::F64vec tvec = tomg ^ tpos;
            fprintf(fp, " %+e", tpot - 0.5 * (tvec * tvec));         // 45
        }
        fprintf(fp, "\n");
    }

    void readHexa(FILE *fp) {
        PS::F64 (*cvt)(PS::U64) = convertU64ToF64;
        PS::S32 s0, s1;
        PS::U64 umass;
        PS::U64 upos[3], uvel[3], uacc[3];
        PS::U64 uuene, uudot;
        PS::U64 ualph, uadot, ualphu, uadotu;
        PS::U64 udens, upres, uvsnd, utemp;
        PS::U64 udivv, urotv, ubswt;
        PS::U64 uksr, ugrdh, uvsmx;
        PS::U64 upot, uenuc;
        PS::U64 ucmps[NR::NumberOfNucleon];

        fscanf(fp, "%d %d %llx", &this->id, &this->istar, &umass);
        fscanf(fp, "%llx %llx %llx", &upos[0], &upos[1], &upos[2]);
        fscanf(fp, "%llx %llx %llx", &uvel[0], &uvel[1], &uvel[2]);
        fscanf(fp, "%llx %llx %llx", &uacc[0], &uacc[1], &uacc[2]);
        fscanf(fp, "%llx %llx", &uuene, &uudot);
        fscanf(fp, "%llx %llx", &ualph, &uadot);
        fscanf(fp, "%llx %llx", &ualphu, &uadotu);
        fscanf(fp, "%llx %llx %llx %llx", &udens, &upres, &uvsnd, &utemp);
        fscanf(fp, "%llx %llx %llx", &udivv, &urotv, &ubswt);
        fscanf(fp, "%llx %llx %llx", &uksr,  &ugrdh, &uvsmx);
        fscanf(fp, "%llx %llx", &upot, &uenuc);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fscanf(fp, "%llx", &ucmps[k]);
        }

        this->mass = cvt(umass);
        this->pos[0] = cvt(upos[0]);
        this->pos[1] = cvt(upos[1]);
        this->pos[2] = cvt(upos[2]);
        this->vel[0] = cvt(uvel[0]);
        this->vel[1] = cvt(uvel[1]);
        this->vel[2] = cvt(uvel[2]);
        this->acc[0] = cvt(uacc[0]);
        this->acc[1] = cvt(uacc[1]);
        this->acc[2] = cvt(uacc[2]);
        this->uene  = cvt(uuene);
        this->udot  = cvt(uudot);
        this->alph  = cvt(ualph);
        this->adot  = cvt(uadot);
        this->alphu = cvt(ualphu);
        this->adotu = cvt(uadotu);
        this->dens  = cvt(udens);
        this->pres  = cvt(upres);
        this->vsnd  = cvt(uvsnd);
        this->temp  = cvt(utemp);
        this->divv  = cvt(udivv);
        this->rotv  = cvt(urotv);
        this->bswt  = cvt(ubswt);
        this->ksr   = cvt(uksr);
        this->grdh  = cvt(ugrdh);
        this->vsmx  = cvt(uvsmx);
        this->pot   = cvt(upot);        
        this->enuc  = cvt(uenuc);        
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            this->cmps[k] = cvt(ucmps[k]);
        }
    }

    void writeHexa(FILE *fp) const {
        PS::U64 (*cvt)(PS::F64) = convertF64ToU64;
        fprintf(fp, "%6d %2d %llx", this->id, this->istar, cvt(this->mass));
        fprintf(fp, " %llx %llx %llx", cvt(this->pos[0]), cvt(this->pos[1]), cvt(this->pos[2]));
        fprintf(fp, " %llx %llx %llx", cvt(this->vel[0]), cvt(this->vel[1]), cvt(this->vel[2]));
        fprintf(fp, " %llx %llx %llx", cvt(this->acc[0]), cvt(this->acc[1]), cvt(this->acc[2]));
        fprintf(fp, " %llx %llx", cvt(this->uene), cvt(this->udot));
        fprintf(fp, " %llx %llx", cvt(this->alph), cvt(this->adot));
        fprintf(fp, " %llx %llx", cvt(this->alphu), cvt(this->adotu));
        fprintf(fp, " %llx %llx", cvt(this->dens), cvt(this->pres));
        fprintf(fp, " %llx %llx", cvt(this->vsnd), cvt(this->temp));
        fprintf(fp, " %llx %llx %llx", cvt(this->divv), cvt(this->rotv), cvt(this->bswt));
        fprintf(fp, " %llx %llx %llx", cvt(this->ksr), cvt(this->grdh), cvt(this->vsmx));
        fprintf(fp, " %llx %llx", cvt(this->pot), cvt(this->enuc));
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fprintf(fp, " %llx", cvt(this->cmps[k]));
        }
        fprintf(fp, "\n");
    }
    
    void referEquationOfState() {
        CalcEquationOfState::getThermodynamicQuantity(this->dens,
                                                      this->uene,
                                                      this->abar,
                                                      this->zbar,
                                                      this->pres,
                                                      this->vsnd,
                                                      this->temp,
                                                      this->cnteos);
    }
    
    void referEquationOfStateDamping1() {
        this->pres = CalcEquationOfState::getPressure(this->dens, this->uene);
        this->vsnd = CalcEquationOfState::getSoundVelocity(this->dens, this->uene);
        this->temp = CalcEquationOfState::getTemperature(this->dens, this->uene);
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

    void calcReleasedNuclearEnergy() {
        if(this->temp > 1e6) {
            this->dnuc  = CalcNRH::getGeneratedEnergy(RP::Timestep,
                                                      this->dens,
                                                      this->temp,
                                                      this->cmps);
            this->enuc += this->dnuc;
        }
    }

    PS::F64 calcTimestep() {
        using namespace CodeUnit;
        PS::F64 tceff   = RP::CoefficientOfTimestep;
        PS::F64 dth = tceff * this->ksr / (this->vsmx * SK::ksrh);
        PS::F64 dtu = tceff * fabs(this->uene / this->udot);
        PS::F64 dtn = tceff * fabs(this->uene / this->dnuc) * RP::Timestep;
        dtu = (this->dens < 1e4 * UnitOfDensity) ? RP::MaximumTimestep : dtu;
        dtn = (this->dens < 1e4 * UnitOfDensity) ? RP::MaximumTimestep : dtn;
        return std::min(std::min(dth, dtu), dtn);
    }

    inline void addAdditionalForceDamping2() {
        this->acc  -= RP::RotationalVelocity ^ (RP::RotationalVelocity ^ this->pos)
            + 2.d * (RP::RotationalVelocity ^ this->vel);
    }

    PS::F64 calcEnergy() {
        return this->mass * (0.5 * this->vel * this->vel + this->uene + 0.5 * this->pot);
    }

    PS::F64 calcEnergyDamping2() {
        PS::F64vec tvec = RP::RotationalVelocity ^ this->pos;
        return this->mass * (0.5 * this->vel * this->vel + this->uene
                             + 0.5 * (this->pot - tvec * tvec));
    }

    PS::F64vec calcMomentum() {
        return this->mass * this->vel;
    }

    void predict(PS::F64 dt) {
        this->pos    = this->pos   +       this->vel   * dt  + 0.5 * this->acc * dt * dt;
        this->vel2   = this->vel   + 0.5 * this->acc   * dt;
        this->vel    = this->vel   +       this->acc   * dt;
        this->uene2  = this->uene  + 0.5 * this->udot  * dt;
        //this->uene   = this->uene  +       this->udot  * dt;
        this->uene   = this->uene  +       this->udot  * dt  +       this->dnuc;
        this->alph2  = this->alph  + 0.5 * this->adot  * dt;
        this->alph   = this->alph  +       this->adot  * dt;
        this->alphu2 = this->alphu + 0.5 * this->adotu * dt;
        this->alphu  = this->alphu +       this->adotu * dt;
    }

    void correct(PS::F64 dt) {
        this->vel   = this->vel2   + 0.5 * this->acc   * dt;
        //this->uene  = this->uene2  + 0.5 * this->udot  * dt;
        this->uene  = this->uene2  + 0.5 * this->udot  * dt  +       this->dnuc;
        this->alph  = this->alph2  + 0.5 * this->adot  * dt;
        this->alphu = this->alphu2 + 0.5 * this->adotu * dt;
    }

    void correctDamping1(PS::F64 dt) {
        this->acc  -= this->vel * 0.05;
        this->vel   = this->vel2   + 0.5 * this->acc   * dt;
        this->uene  = this->uene2  + 0.5 * this->udot  * dt;
        this->alph  = this->alph2  + 0.5 * this->adot  * dt;
        this->alphu = this->alphu2 + 0.5 * this->adotu * dt;
        this->uene  = CalcEquationOfState::getEnergy(this->dens, this->uene);
    }    

    void correctDamping2(PS::F64 dt) {
        this->acc  -= this->vel * RP::ReductionTimeInv;
        this->vel   = this->vel2   + 0.5 * this->acc  * dt;
        this->uene  = this->uene2  + 0.5 * this->udot * dt;
        this->alph  = this->alph2  + 0.5 * this->adot * dt;
        this->alphu = this->alphu2 + 0.5 * this->adotu * dt;
    }

    void correctDamping3(PS::F64 dt) {
        this->vel   = this->vel2   + 0.5 * this->acc   * dt;
        this->uene  = this->uene2  + 0.5 * this->udot  * dt;
        this->alph  = this->alph2  + 0.5 * this->adot  * dt;
        this->alphu = this->alphu2 + 0.5 * this->adotu * dt;
        PS::F64 unow = this->uene;
        PS::F64 umin = CalcEquationOfState::getEnergyMin(this->dens, this->abar, this->zbar);
        this->uene   = (unow < umin) ? unow : ((unow - umin) * exp(-0.1 * dt) + umin);
    }

    void calcAlphaDot() {
        PS::F64 src    = std::max((- divv * (RP::AlphaMaximum - this->alph)), 0.);
        PS::F64 tauinv = (0.25 * SK::ksrh * this->vsnd) / this->ksr;
        this->adot     = - (this->alph - RP::AlphaMinimum) * tauinv + src;

        PS::F64 uene   = std::max(this->uene - CalcEquationOfState::getEnergyMin(this->dens, this->abar, this->zbar), 0.);
        PS::F64 srcu   = this->ksr * SK::ksrhinv * fabs(this->diffu)
            / sqrt(uene + RP::EpsilonOfInternalEnergy)
            * std::max(RP::AlphuMaximum - this->alphu, 0.);
        this->adotu    = - (this->alphu - RP::AlphuMinimum) * tauinv + src;
    }
};

namespace RunParameter {    

    PS::S64 readAscii(const char * const filename) {
        using namespace CodeUnit;

        FILE *fp = fopen(filename, "r");
        if(fp == NULL) {
            fprintf(stderr, "Not found header file %s!\n", filename);
            PS::Abort();
        }
        fscanf(fp, "%lf%lf%lf%lf", &Time, &TimeEnd, &TimestepAscii, &TimestepHexa);
        fscanf(fp, "%lld%lld", &FlagDamping, &FlagNuclear);
        fscanf(fp, "%lf%lf", &AlphaMaximum, &AlphaMinimum);
        fscanf(fp, "%lf%lf", &AlphuMaximum, &AlphuMinimum);
        fscanf(fp, "%lf", &CoefficientOfTimestep);
        fscanf(fp, "%lld", &NumberOfParticle);
        fclose(fp);

        Time          *= UnitOfTimeInv;
        TimeEnd       *= UnitOfTimeInv;
        TimestepAscii *= UnitOfTimeInv;
        TimestepHexa  *= UnitOfTimeInv;
    }


    template <class Tdinfo,
              class Tsph,
              class Tmsls>
    void readHexa(FILE *fp,                         
                  Tdinfo & dinfo,
                  Tsph & sph,
                  Tmsls & msls) {
        PS::F64 (*cvt)(PS::U64) = convertU64ToF64;        
        PS::U64 utime, utend, udtime, udta, udth;
        PS::U64 ualphamax, ualphamin, ualphumax, ualphumin;
        PS::U64 uksrmax, uepsu, uenuc;
        PS::U64 urtime, urtimeinv, utcoeff;
        PS::U64 urv[3];
        PS::S64 nstep, nptcl, nmsls, fdamp, fnuc;
        fscanf(fp, "%llx %llx %llx %llx %llx %lld", &utime, &utend, &udtime, &udta, &udth, &nstep);
        fscanf(fp, "%llx %llx %llx %llx", &ualphamax, &ualphamin, &ualphumax, &ualphumin);
        fscanf(fp, "%llx %llx", &uksrmax, &uepsu);
        fscanf(fp, "%llx %llx", &urtime, &urtimeinv);
        fscanf(fp, "%lld %lld %llx", &nptcl, &nmsls, &utcoeff);
        fscanf(fp, "%llx %llx %llx", &urv[0], &urv[1], &urv[2]);
        fscanf(fp, "%llx %lld %lld", &uenuc, &fdamp, &fnuc);
        Time          = cvt(utime);
        TimeEnd       = cvt(utend);
        Timestep      = cvt(udtime);
        TimestepAscii = cvt(udta);
        TimestepHexa  = cvt(udth);
        NumberOfStep  = nstep;
        AlphaMaximum  = cvt(ualphamax);
        AlphaMinimum  = cvt(ualphamin);
        AlphuMaximum  = cvt(ualphumax);
        AlphuMinimum  = cvt(ualphumin);
        KernelSupportRadiusMaximum = cvt(uksrmax);
        EpsilonOfInternalEnergy    = cvt(uepsu);
        ReductionTime    = cvt(urtime);
        ReductionTimeInv = cvt(urtimeinv);
        sph.setNumberOfParticleLocal(nptcl);
        msls.setNumberOfParticleLocal(nmsls);
        CoefficientOfTimestep = cvt(utcoeff);
        RotationalVelocity[0] = cvt(urv[0]);
        RotationalVelocity[1] = cvt(urv[1]);
        RotationalVelocity[2] = cvt(urv[2]);
        NuclearEnergyTotal    = cvt(uenuc);
        FlagDamping           = fdamp;
        FlagNuclear           = fnuc;

        for(PS::S32 i = 0; i < PS::Comm::getNumberOfProc(); i++) {
            PS::F64ort domain;
            PS::U64    u0, u1, u2, u3, u4, u5;
            fscanf(fp, "%llx %llx %llx %llx %llx %llx\n", &u0, &u1, &u2, &u3, &u4, &u5);
            domain.low_[0]  = cvt(u0);
            domain.low_[1]  = cvt(u1);
            domain.low_[2]  = cvt(u2);
            domain.high_[0] = cvt(u3);
            domain.high_[1] = cvt(u4);
            domain.high_[2] = cvt(u5);
            dinfo.setPosDomain(i, domain);
        }
    }

    template <class Tdinfo,
              class Tsph,
              class Tmsls>
    void writeHexa(FILE *fp,                          
                   Tdinfo & dinfo,
                   Tsph & sph,
                   Tmsls & msls) {
        PS::U64 (*cvt)(PS::F64) = convertF64ToU64;
        fprintf(fp, "%llx %llx %llx %llx %llx %lld\n",
                cvt(Time), cvt(TimeEnd), cvt(Timestep),
                cvt(TimestepAscii), cvt(TimestepHexa), NumberOfStep);
        fprintf(fp, "%llx %llx %llx %llx\n",
                cvt(AlphaMaximum), cvt(AlphaMinimum), cvt(AlphuMaximum), cvt(AlphuMinimum));
        fprintf(fp, "%llx %llx %llx %llx\n",
                cvt(KernelSupportRadiusMaximum), cvt(EpsilonOfInternalEnergy),
                cvt(ReductionTime), cvt(ReductionTimeInv));
        fprintf(fp, "%lld %lld %llx\n", sph.getNumberOfParticleLocal(),
                msls.getNumberOfParticleLocal(), cvt(CoefficientOfTimestep));
        fprintf(fp, "%llx %llx %llx\n", cvt(RotationalVelocity[0]),
                cvt(RotationalVelocity[1]), cvt(RotationalVelocity[2]));
        fprintf(fp, "%llx %lld %lld\n", cvt(NuclearEnergyTotal), FlagDamping, FlagNuclear);

        PS::S32 nproc = PS::Comm::getNumberOfProc();
        for(PS::S32 i = 0; i < nproc; i++) {
            PS::F64ort domain = dinfo.getPosDomain(i);
            fprintf(fp, "%llx %llx %llx %llx %llx %llx\n",
                    cvt(domain.low_[0]),  cvt(domain.low_[1]),  cvt(domain.low_[2]),
                    cvt(domain.high_[0]), cvt(domain.high_[1]), cvt(domain.high_[2]));
        }
    }

    void outputRunParameter(char **argv) {
        if(PS::Comm::getRank() != 0) {
            return;
        }
        FILE * fp = FilePointerForLog;
        fprintf(fp, "# # of process %8d\n", PS::Comm::getNumberOfProc());
        fprintf(fp, "# # of thread  %8d\n", PS::Comm::getNumberOfThread());
        if(atoi(argv[1]) == 0) {
            fprintf(fp, "# Header Data:  %s\n", argv[3]);
        } else {
            fprintf(fp, "# Restart file: %s\n", argv[2]);
        }
        fprintf(fp, "# AlphaMax:   %+e\n", AlphaMaximum);
        fprintf(fp, "# AlphaMin:   %+e\n", AlphaMinimum);
        fprintf(fp, "# AlphuMax:   %+e\n", AlphuMaximum);
        fprintf(fp, "# AlphuMin:   %+e\n", AlphuMinimum);
        fprintf(fp, "# C for Time: %+e\n", CoefficientOfTimestep);
        fprintf(fp, "# KernelMax:  %+e\n", KernelSupportRadiusMaximum * CodeUnit::UnitOfLength);
        fprintf(fp, "# EpsilonU:   %+e\n", EpsilonOfInternalEnergy * CodeUnit::UnitOfEnergy);
        fprintf(fp, "# Damping: %d\n", FlagDamping);
        fprintf(fp, "# Nuclear: %d\n", FlagNuclear);        
        fprintf(fp, "# # of Step %8d\n", NumberOfStep);
        fflush(fp);
    }

    template <class Tsph>
    PS::F64 setEpsilonOfInternalEnergy(Tsph & sph) {
        PS::F64 uminloc = std::numeric_limits<double>::max();
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].uene < uminloc) {
                uminloc = sph[i].uene;
            }
        }
        PS::F64 uminglb = PS::Comm::getMinValue(uminloc);
        return 1e-4 * uminglb;
    }
    
}

template <class Tdinfo,
          class Tsph,
          class Tmsls>
void readHexa(char * const filename,
              Tdinfo & dinfo,
              Tsph & sph,
              Tmsls & msls) {
    FILE *fp = fopen(filename, "r");
    RP::readHexa(fp, dinfo, sph, msls);
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].readHexa(fp);
    }
    for(PS::S32 i = 0; i < msls.getNumberOfParticleLocal(); i++) {
        msls[i].readHexa(fp);
    }
    fclose(fp);
    RP::NumberOfParticle = PS::Comm::getSum(sph.getNumberOfParticleLocal());
}

template <class Tdinfo,
          class Tsph,
          class Tmsls>
void writeHexa(char * const filename,
               Tdinfo & dinfo,
               Tsph & sph,
               Tmsls & msls) {

    FILE *fp = fopen(filename, "w");
    RP::writeHexa(fp, dinfo, sph, msls);
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].writeHexa(fp);
    }
    for(PS::S32 i = 0; i < msls.getNumberOfParticleLocal(); i++) {
        msls[i].writeHexa(fp);
    }
    fclose(fp);
}

template <class Tdinfo,
          class Tsph,
          class Tmsls>
void outputData(Tdinfo & dinfo,
                Tsph & sph,
                Tmsls & msls) {
    if(RP::Time - (PS::S64)(RP::Time / RP::TimestepAscii) * RP::TimestepAscii == 0.) {
        char filename[64];
        sprintf(filename, "snap/sph_t%04d.dat", (PS::S32)RP::Time);
        sph.writeParticleAscii(filename);
        if(RP::FlagDamping == 2) {
            sprintf(filename, "snap/msls_t%04d.dat", (PS::S32)RP::Time);
            msls.writeParticleAscii(filename);
        }
    }
    if(RP::Time - (PS::S64)(RP::Time / RP::TimestepHexa) * RP::TimestepHexa == 0.) {
        char filename[64];
        sprintf(filename, "snap/t%04d_p%06d.hexa", (PS::S32)RP::Time, PS::Comm::getRank());
        writeHexa(filename, dinfo, sph, msls);
    }
    PS::F64 etot = calcEnergy(sph);
    PS::F64 enuc = calcReleasedNuclearEnergyTotal(sph);
    if(PS::Comm::getRank() == 0) {
        using namespace CodeUnit;
        fprintf(RP::FilePointerForLog, "time: %8d %16.10f %+e %+e %+e %+e %+e\n",
                RP::NumberOfStep, RP::Time * UnitOfTime, RP::Timestep * UnitOfTime,
                (etot - enuc) * UnitOfEnergy * UnitOfMass, etot * UnitOfEnergy * UnitOfMass,
                enuc * UnitOfEnergy * UnitOfMass, WT::getTimeTotal());
        WT::dump(RP::Time, RP::FilePointerForTime);
        fflush(RP::FilePointerForLog);
    }
}

template <class Tsph>
void predict(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].predict(RP::Timestep);
    }
}

template <class Tsph>
void correct(Tsph & sph) {
    if(RP::FlagDamping == 0) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].correct(RP::Timestep);
        }
    } else if(RP::FlagDamping == 1) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].correctDamping1(RP::Timestep);
        }
    } else if(RP::FlagDamping == 2) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].correctDamping2(RP::Timestep);
        }
    } else if(RP::FlagDamping == 3) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].correctDamping3(RP::Timestep);
        }
    } else {
        fprintf(stderr, "Not supported damping mode %d!\n", RP::FlagDamping);
        PS::Abort();
    }
}

void initializeSimulation() {
    using namespace RunParameter;
    //Time              = 0.;
    //TimeEnd           = 0.;
    MaximumTimestep   = 1. / 64.;
    MinimumTimestep   = 1e-16;
    Timestep          = RP::MaximumTimestep;
    //TimeAscii         = 0.;
    //TimestepAscii     = 0.;
    //TimeHexa          = 0.;
    //TimestepHexa      = 0.;
    NumberOfStep      = 0;
    //NumberOfAscii     = 0;
    //NumberOfHexa      = 0;
    ComputationalBox.high_ = 0;
    ComputationalBox.low_  = 0;
    NumberOfDimension = 3;
    //KernelType        = 1;
    //CoefficientOfTimestep      = 0.1;
    GravitationalSoftening     = 3e6 * CodeUnit::UnitOfLengthInv;
    KernelSupportRadiusMaximum = 0.;
    EpsilonOfInternalEnergy    = 0.;    
    //NumberOfParticle   = 0;
    //AdiabaticIndex;
    //RotationalVelocity; // omega
    NuclearEnergyTotal = 0.; // enuc
    FlagGravity = 1;
    //FlagDamping = 0;
    //FlagNuclear = 0;
    if(PS::Comm::getRank() == 0) {
        RP::FilePointerForLog  = fopen("snap/time.log", "w");
        RP::FilePointerForTime = fopen("snap/prof.log", "w");
    }
}

template <class Tdinfo,
          class Tsph,
          class Tmsls,
          class Tdensity,
          class Thydro,
          class Tgravity>
void startSimulation(char **argv,
                     Tdinfo & dinfo,
                     Tsph & sph,
                     Tmsls & msls,
                     Tdensity & density,
                     Thydro & hydro,
                     Tgravity & gravity) {
    char headname[256];
    sprintf(headname, "%s.head", argv[3]);
    RP::readAscii(headname);
    if(atoi(argv[2]) == 0) {
        char filename[256];
        sprintf(filename, "%s.data", argv[3]);
        sph.readParticleAscii(filename);    
    } else {
        sph.readParticleAscii(argv[3], "%s_p%06d_i%06d.data");    
    }
    if(RP::FlagDamping == 2) {
        generateMassLessParticle(msls, sph);
    }
    ND::setDimension(RP::NumberOfDimension);
    SK::setKernel(RP::KernelType, RP::NumberOfDimension);
    RP::KernelSupportRadiusMaximum = calcSystemSize(sph);
    RP::EpsilonOfInternalEnergy    = RP::setEpsilonOfInternalEnergy(sph);
    RP::outputRunParameter(argv);
    PS::MT::init_genrand(0);
    dinfo.decomposeDomainAll(sph);
    sph.exchangeParticle(dinfo);
    if(RP::FlagDamping == 2) {
        msls.exchangeParticle(dinfo);
        calcRotationalVelocity(sph);
    }
    calcSPHKernel(dinfo, sph, msls, density, hydro, gravity);
}

template <class Tdinfo,
          class Tsph,
          class Tmsls>
void restartSimulation(char **argv,
                       Tdinfo & dinfo,
                       Tsph & sph,
                       Tmsls & msls) {
    char filename[64];
    sprintf(filename, "%s_p%06d.hexa", argv[2], PS::Comm::getRank());
    readHexa(filename, dinfo, sph, msls);
    ND::setDimension(RP::NumberOfDimension);
    SK::setKernel(RP::KernelType, RP::NumberOfDimension);
    RP::outputRunParameter(argv);
    return;
}

template <class Tdinfo,
          class Tsph,
          class Tmsls,
          class Tdensity,
          class Thydro,
          class Tgravity>
void loopSimulation(Tdinfo & dinfo,
                    Tsph & sph,
                    Tmsls & msls,
                    Tdensity & density,
                    Thydro & hydro,
                    Tgravity & gravity) {

    WT::clear();
    while(RP::Time < RP::TimeEnd) {
        RP::Timestep = calcTimestep(sph);
        WT::reduceInterProcess();
        outputData(dinfo, sph, msls);
        WT::clear();
        if(RP::FlagDamping == 2) {
            reduceSeparation(sph, msls);
        }
        WT::start();
        if(RP::FlagNuclear == 1) {
            calcReleasedNuclearEnergy(sph);
        }
        WT::accumulateCalcNuclearReaction();
        WT::start();
        predict(sph);
        WT::accumulateOthers();
        WT::start();
        if(RP::NumberOfStep % 4 == 0) {
            PS::MT::init_genrand(0);
            dinfo.decomposeDomainAll(sph);
        }
        WT::accumulateDecomposeDomain();
        WT::start();
        sph.exchangeParticle(dinfo);
        WT::accumulateExchangeParticle();
        calcSPHKernel(dinfo, sph, msls, density, hydro, gravity);
        WT::start();
        correct(sph);
        WT::accumulateOthers();
        RP::Time += RP::Timestep;
        RP::NumberOfStep++;
    }
}

template <class Tdinfo,
          class Tsph,
          class Tmsls>
void finalizeSimulation(Tdinfo & dinfo,
                        Tsph & sph,
                        Tmsls & msls) {
    outputData(dinfo, sph, msls);
    if(PS::Comm::getRank() == 0) {
        fclose(RP::FilePointerForLog);
        fclose(RP::FilePointerForTime);
    }
    PS::F64    mc;
    PS::F64vec xc;
    PS::F64vec vc;
    calcCenterOfMass(sph, mc, xc, vc);
    PS::S32 nloc = sph.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        sph[i].pos -= xc;
        sph[i].vel -= vc;
    }
    sph.writeParticleAscii("snap/final.dat");
}

typedef HelmholtzGas GeneralSPH;
