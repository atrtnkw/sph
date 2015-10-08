#pragma once

class Density;
class DensityEPI;
class DensityEPJ;
class Volume;
class Pressure;
class Auxiliary;
class Hydro;
class Gravity;

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
    PS::F64    alphu;
    PS::F64    alphu2;    
    PS::F64    adotu;
    PS::F64    diffu;
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

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {        
        this->pos = pos_new;
    }


    void copyFromForce(const Density & density);
    void copyFromForce(const Volume & volume);
    void copyFromForce(const Pressure & pressure);
    void copyFromForce(const Auxiliary & auxiliary);
    void copyFromForce(const Hydro & hydro);
    void copyFromForce(const Gravity & gravity);

    virtual void readAscii(FILE * fp) {};
    virtual void writeAscii(FILE *fp) const {};
    virtual void referEquationOfState() {};
    virtual PS::F64 calcTimestep() {};
    virtual PS::F64 calcEnergy() {};
    virtual PS::F64 calcEnergyDamping2() {};
    virtual void predict(PS::F64 dt) {};
    virtual void correct(PS::F64 dt) {};
    virtual void calcAlphaDot() {};
    virtual void addAdditionalForce() {};
    virtual void addAdditionalForceDamping2() {};

    void calcBalsaraSwitch() {
        this->bswt = fabs(this->divv)
            / (fabs(this->divv) + fabs(this->rotv) + 1e-4 * this->vsnd * SK::ksrh / this->ksr);
    }

    void copyAcceleration() {
        this->acc = this->acch;
    }

    void sumAcceleration() {
        this->acc = this->acch + this->accg1 + this->accg2;
    }

    //void readHexa(FILE *fp);
    //void writeHexa(FILE *fp) const;

};

template <class Tsph>
void referEquationOfState(Tsph & sph);
template <class Tsph>
void calcBalsaraSwitch(Tsph & sph);
template <class Tsph>
void calcAlphaDot(Tsph & sph);
template <class Tsph>
PS::F64 calcEnergy(Tsph & sph);
template <class Tsph>
PS::F64vec calcMomentum(Tsph & sph);
