#pragma once

#include "hdr_heos.hpp"
#include "hdr_nuc.hpp"

class Header {
public:
    PS::S32 nptcl;
    PS::F64vec cpos;
    PS::S64 msls_nx;
    PS::F64 msls_xmin;
    PS::F64 msls_xmax;

    Header() {
        nptcl     = 0;
        cpos      = 0.;
        msls_nx   = 0;
        msls_xmin = 0.;
        msls_xmax = 0.;
    }

    void initialize(char * filename) {
        FILE *fp = fopen(filename, "r");
        fscanf(fp, "%d",  &this->nptcl);
        fscanf(fp, "%lf%lf%lf", &this->cpos[0], &this->cpos[1], &this->cpos[2]);
        fscanf(fp, "%d",  &this->msls_nx);
        fscanf(fp, "%lf", &this->msls_xmin);
        fscanf(fp, "%lf", &this->msls_xmax);
        fclose(fp);
        this->cpos      *= CodeUnit::UnitOfLengthInv;
        this->msls_xmin *= CodeUnit::UnitOfLengthInv;
        this->msls_xmax *= CodeUnit::UnitOfLengthInv;
    }

    PS::S32 readAscii(FILE *fp) {
        return this->nptcl;
    }
};

class AnalysisGas : public SPH {
public:
    PS::S64 istar;
    PS::F64 temp;
    PS::F64 cmps[NR::NumberOfNucleon];
    PS::S64 cnteos;    
    PS::F64 dnuc;
    PS::F64 enuc;
    bool    fnse;

    void readAscii(FILE *fp) {
        using namespace CodeUnit;

        fscanf(fp, "%d%d%lf",   &this->id,     &this->istar,  &this->mass);
        fscanf(fp, "%lf%lf%lf", &this->pos[0], &this->pos[1], &this->pos[2]);
        fscanf(fp, "%lf%lf%lf", &this->vel[0], &this->vel[1], &this->vel[2]);
        fscanf(fp, "%lf%lf%lf", &this->acc[0], &this->acc[1], &this->acc[2]);
        fscanf(fp, "%lf%lf%lf", &this->uene,   &this->alph,   &this->alphu);
        fscanf(fp, "%lf%lf%d",  &this->dens,   &this->ksr,    &this->np);
        fscanf(fp, "%lf%lf%lf", &this->vsnd,   &this->pres,   &this->temp);
        fscanf(fp, "%lf%lf%lf", &this->divv,   &this->rotv,   &this->bswt);
        fscanf(fp, "%lf%lf%lf", &this->pot,    &this->abar,   &this->zbar);
        fscanf(fp, "%lf",       &this->enuc);
        fscanf(fp, "%lf%lf%lf", &this->vsmx,   &this->udot,   &this->dnuc);
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) {
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        //fscanf(fp, "%d", &this->fnse);

        this->mass *= UnitOfMassInv;
        this->pos  *= UnitOfLengthInv;
        this->vel  *= UnitOfVelocityInv;
        this->acc  /= UnitOfAcceleration;
        this->uene *= UnitOfEnergyInv;
        this->dens *= UnitOfDensityInv;
        this->ksr  *= UnitOfLengthInv;
        this->vsnd *= UnitOfVelocityInv;
        this->pres *= UnitOfPressureInv;
        this->divv *= UnitOfTime;
        this->rotv *= UnitOfTime;
        this->pot  *= UnitOfEnergyInv;
        this->enuc *= UnitOfEnergyInv;
        this->vsmx *= UnitOfVelocityInv;
        this->udot *= UnitOfEnergyInv * UnitOfTime;
        this->dnuc *= UnitOfEnergyInv;
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
        //fprintf(fp, " %d", this->fnse);
        fprintf(fp, "\n");
    }

};

typedef AnalysisGas GeneralSPH;
