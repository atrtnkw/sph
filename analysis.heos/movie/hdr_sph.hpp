#pragma once

const PS::S64 NumberOfNucleon = 13;

class Quantity{
public:
    PS::F64 dens;
    PS::F64 temp;
    PS::F64 abar;
    void clear(){
        dens = 0.0;
        temp = 0.0;
        abar = 0.0;
    }
};

class Header {
public:
    PS::S32 nptcl;
    Header() {
        nptcl = 0;
    }
    PS::S32 readAscii(FILE *fp) {
        fscanf(fp, "%d", &this->nptcl);
        return this->nptcl;
    }
};

class SPH{
public:
    PS::S64    id;
    PS::S64    istar;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    uene;
    PS::F64    alph;
    PS::F64    alphu;
    PS::F64    ksr;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    pres;
    PS::F64    temp;
    PS::F64    divv;
    PS::F64    rotv;
    PS::F64    bswt;
    PS::F64    grdh;
    PS::S64    np;
    PS::F64    pot;
    PS::F64    abar;
    PS::F64    zbar;
    PS::F64    cmps[NumberOfNucleon];
    PS::F64    vsmx;
    PS::F64    udot;
    PS::F64    dnuc;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {        
        this->pos = pos_new;
    }

    void readAscii(FILE *fp) {
        fscanf(fp, "%lld%lld%lf", &this->id, &this->istar, &this->mass);
        fscanf(fp, "%lf%lf%lf", &this->pos[0], &this->pos[1], &this->pos[2]);
        fscanf(fp, "%lf%lf%lf", &this->vel[0], &this->vel[1], &this->vel[2]);
        fscanf(fp, "%lf%lf%lf", &this->acc[0], &this->acc[1], &this->acc[2]);
        fscanf(fp, "%lf%lf%lf%lf", &this->uene, &this->alph, &this->alphu, &this->ksr);
        fscanf(fp, "%lf%lf%lf%lf", &this->dens, &this->vsnd, &this->pres, &this->temp);
        fscanf(fp, "%lf%lf%lf", &this->divv, &this->rotv, &this->bswt);
        fscanf(fp, "%lf%lld%lf", &this->grdh, &this->np, &this->pot);
        fscanf(fp, "%lf%lf", &this->abar, &this->zbar);
        for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        fscanf(fp, "%lf", &this->vsmx);
        fscanf(fp, "%lf", &this->udot);
        fscanf(fp, "%lf", &this->dnuc);
    }

    void writeAscii(FILE *fp) {
        fprintf(fp, "%+e %+e %+e\n", this->pos[0], this->pos[1], this->pos[2]);
    }

    static inline PS::F64 calcVolumeInverse(const PS::F64 hi);
    static inline PS::F64 calcPowerOfDimInverse(PS::F64 mass,
                                                PS::F64 dens);
    static inline v4df calcVolumeInverse(const v4df hi);
    static inline v8sf calcVolumeInverse(const v8sf hi);

};

inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi * hi * hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return pow(mass / dens, 1.d / 3.d);
}
inline v4df SPH::calcVolumeInverse(const v4df hi) {return hi * hi * hi;}
inline v8sf SPH::calcVolumeInverse(const v8sf hi) {return hi * hi * hi;}
