#pragma once

class Density;
class DensityEPI;
class DensityEPJ;
class Volume;
class Auxiliary;
class Hydro;

class SPH {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec vel2;
    PS::F64vec acc;
    PS::F64    uene;
    PS::F64    uene2;
    PS::F64    udot;
    PS::F64    alph;
    PS::F64    alph2;
    PS::F64    adot;
    PS::F64    vol;
    PS::F64    ksr;
    PS::F64    rs;
    PS::F64    dens;
    PS::F64    pres;
    PS::F64    vsnd;
    PS::F64    divv;
    PS::F64    rotv;
    PS::F64    bswt;
    PS::S64    np;
    PS::F64    grdh;
    PS::F64    vsmx;
    PS::F64    pot;
    PS::F64vec acch;
    PS::F64vec accg1;
    PS::F64vec accg2;
    PS::F64    eta;
    PS::F64    abar;
    PS::F64    zbar;
    PS::F64mat ctau;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {        
        this->pos = pos_new;
    }


    void copyFromForce(const Density & density);
    void copyFromForce(const Volume & volume);
    void copyFromForce(const Auxiliary & auxiliary);
    void copyFromForce(const Hydro & hydro);

    virtual void readAscii(FILE * fp) {};
    virtual void writeAscii(FILE *fp) const {};
    virtual void referEquationOfState() {};

    void calcBalsaraSwitch() {
        this->bswt = fabs(this->divv)
            / (fabs(this->divv) + fabs(this->rotv) + 1e-4 * this->vsnd * SK::ksrh / this->ksr);
    }

    //void calcAlphaDot();
    //PS::F64 calcTimeStep();
    //PS::F64 calcEnergy();
    //void predict(PS::F64 dt) {};
    //void correct(PS::F64 dt) {};
    //void calcAbarZbar();
    //inline void addAdditionalForce();
    //void readHexa(FILE *fp);
    //void writeHexa(FILE *fp) const;

    //static inline PS::F64 calcVolumeInverse(const PS::F64 hi);
    //static inline PS::F64 calcPowerOfDimInverse(PS::F64 mass,
    //PS::F64 dens);
    //static inline v4df calcVolumeInverse(const v4df hi);
    //static inline v8sf calcVolumeInverse(const v8sf hi);

};
