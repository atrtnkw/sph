#pragma once

#include "heos/WDEOS.hpp"

namespace CodeUnit {
    PS::F64 SolarRadius        = 6.9599e10; // cm
    PS::F64 SolarMass          = 1.9891e33; // g
    PS::F64 GravityConstant    = 6.6738480e-8;  // [cm^3 g^-1 g^-2]

    PS::F64 UnitOfLength       = SolarRadius * 1.0e-3d;
    PS::F64 UnitOfMass         = SolarMass   * 1.0e-6d;
    PS::F64 UnitOfTime         = 1.0d;

    PS::F64 UnitOfVelocity     = UnitOfLength / UnitOfTime;
    PS::F64 UnitOfDensity      = UnitOfMass / (UnitOfLength * UnitOfLength * UnitOfLength);
    PS::F64 UnitOfEnergy       = UnitOfVelocity * UnitOfVelocity;
    PS::F64 UnitOfPressure     = UnitOfMass * UnitOfLength
        / (UnitOfTime * UnitOfTime * UnitOfLength * UnitOfLength);
    PS::F64 UnitOfAcceleration = UnitOfVelocity / UnitOfTime;

    PS::F64 UnitOfLengthInv    = 1.d / UnitOfLength;
    PS::F64 UnitOfMassInv      = 1.d / UnitOfMass;
    PS::F64 UnitOfTimeInv      = 1.d / UnitOfTime;
    PS::F64 UnitOfVelocityInv  = 1.d / UnitOfVelocity;
    PS::F64 UnitOfDensityInv   = 1.d / UnitOfDensity;
    PS::F64 UnitOfEnergyInv    = 1.d / UnitOfEnergy;
    PS::F64 UnitOfPressureInv  = 1.d / UnitOfPressure;

    // different from OTOO code in time unit !!
    PS::F64 grav               = GravityConstant /
        ((UnitOfLength * UnitOfLength * UnitOfLength) * UnitOfMassInv
         * (UnitOfTimeInv * UnitOfTimeInv));
}

class CalcEquationOfState {
private:
    OTOO::EOS * eos_;
    CalcEquationOfState() {
        using namespace OTOO;
        using namespace CodeUnit;
#ifdef WD_DAMPING
        eos_ = new WDEOSforDumping(UnitOfDensity, UnitOfEnergy, UnitOfVelocity,
                                   UnitOfPressure, 1.0e6);
#else
        eos_ = new WDEOS_D_E(UnitOfDensity, UnitOfEnergy, UnitOfVelocity,
                             UnitOfPressure, 1.0e6);
#endif
    };
    ~CalcEquationOfState() {};
    CalcEquationOfState(const CalcEquationOfState & c);
    CalcEquationOfState & operator=(const CalcEquationOfState & c);
    static CalcEquationOfState & getInstance() {
        static CalcEquationOfState inst;
        return inst;
    }
public:
    static PS::F64 getPressure(PS::F64 density,
                               PS::F64 energy) {
        return getInstance().eos_->GetP(density, energy);
    }
    static PS::F64 getTemperature(PS::F64 density,
                                  PS::F64 energy) {
        return getInstance().eos_->GetT(density, energy);
    }
    static PS::F64 getSoundVelocity(PS::F64 density,
                                    PS::F64 energy) {
        return getInstance().eos_->GetS(density, energy);
    }
    static PS::F64 getEnergy(PS::F64 density,
                             PS::F64 energy) {
        return getInstance().eos_->GetE(density, energy);
    }
};

class Density{
public:
    PS::F64 dens;
    PS::F64 grdh;
    PS::F64 ksr;
    PS::F64 rotv;
    PS::F64 divv;
    PS::S32 np;
    bool    itr;
    void clear(){
        dens = 0.0;
        grdh = 0.0;
        ksr  = 0.0;
        rotv = 0.0;
        divv = 0.0;
        np   = 0;
        itr  = false;
    }
};

class Derivative{
public:
    PS::F64vec acc;
    PS::F64    udot;
    PS::F64    vsmx;
    void clear(){
        acc  = 0.0;
        udot = 0.0;
        vsmx = 0.0;
    }
};

class Gravity{
public:
    PS::F64vec acc;
    PS::F64    pot;
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
};

class Header {
public:
    PS::F64    time;
    PS::F64    tend;
    PS::F64    dtsp;
    PS::F64ort cbox;
    PS::F64    alphamax;
    PS::F64    alphamin;
    PS::F64    tceff;
    PS::F64    eps;
    PS::S64    nptcl;

    Header() {
        time       = 0.0;
        tend       = 0.0;
        dtsp       = 0.0;
        cbox.low_  = 0.0;
        cbox.high_ = 0.0;
        alphamax   = 0.0;
        alphamin   = 0.0;
        tceff      = 0.0;
        eps        = 0.0;
        nptcl      = 0;
    }

    PS::S32 readAscii(FILE *fp) {
        using namespace CodeUnit;

        fscanf(fp, "%lf%lf%lf", &time, &tend, &dtsp);
        fscanf(fp, "%lf%lf%lf", &alphamax, &alphamin, &tceff);
        fscanf(fp, "%lf", &eps);
        fscanf(fp, "%d", &nptcl);

        time *= UnitOfTimeInv;
        tend *= UnitOfTimeInv;
        dtsp *= UnitOfTimeInv;
        eps  *= UnitOfLengthInv;

        return nptcl;
    }
};

class SPH{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    uene;
    PS::F64    udot;
    PS::F64    alph;
    PS::F64    alph2;
    PS::F64    adot;
    PS::F64    dens;
    PS::F64    pres;
    PS::F64    vsnd;
    PS::F64    divv;
    PS::F64    rotv;
    PS::F64    bswt;
    PS::S64    np;
    PS::F64vec vel2;
    PS::F64    uene2;
    PS::F64    ksr;
    PS::F64    rs;
    PS::F64    grdh;
    PS::F64    vsmx;
    PS::F64    pot;
    PS::F64vec accg;
    PS::F64    temp;
    static PS::F64ort cbox;
    static PS::F64    cinv;
    static PS::F64    alphamax, alphamin;
    static PS::F64    tceff;
    static PS::F64    eps;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {        
        this->pos = pos_new;
    }

    void copyFromForce(const Density & density){
        this->dens = density.dens;
        this->grdh = density.grdh;
        this->np   = density.np;
        this->ksr  = density.ksr;
        this->rotv = density.rotv;
        this->divv = density.divv;
    }

    void copyFromForce(const Derivative & derivative){
        this->acc  = derivative.acc;
        this->udot = derivative.udot;
        this->vsmx = derivative.vsmx;
    }

    void copyFromForce(const Gravity & gravity) {
        this->acc  += CodeUnit::grav * gravity.acc;
        this->pot   = CodeUnit::grav * (gravity.pot + this->mass / this->eps);
        this->accg  = CodeUnit::grav * gravity.acc;
    }

    void readAscii(FILE *fp) {        
        using namespace CodeUnit;

        fscanf(fp, "%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->mass,
               &this->pos[0], &this->pos[1], &this->pos[2],
               &this->vel[0], &this->vel[1], &this->vel[2],
               &this->uene,   &this->alph,   &this->ksr);

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
        
        fprintf(fp, "%6d %+e", this->id, tmass);
        fprintf(fp, " %+e %+e %+e", tpos[0], tpos[1], tpos[2]);
        fprintf(fp, " %+e %+e %+e", tvel[0], tvel[1], tvel[2]);
        fprintf(fp, " %+e %+e %+e", tacc[0], tacc[1], tacc[2]);
        fprintf(fp, " %+e %+e %+e", tuene, this->alph, tksr);
        fprintf(fp, " %+e %+e %+e %+e", tdens, tvsnd, tpres, this->temp);
        fprintf(fp, " %+e %+e %+e", tdivv, trotv, this->bswt);
        fprintf(fp, " %+e %6d %+e", this->grdh, this->np, tpot);
        fprintf(fp, "\n");

    }

    void referEquationOfState() {
        this->pres = CalcEquationOfState::getPressure(this->dens, this->uene);
        this->vsnd = CalcEquationOfState::getSoundVelocity(this->dens, this->uene);
        this->temp = CalcEquationOfState::getTemperature(this->dens, this->uene);
    }

    void calcBalsaraSwitch() {
        this->bswt = fabs(this->grdh * this->divv)
            / (fabs(this->grdh * this->divv) + fabs(this->grdh * this->rotv)
               + 1e-4 * this->vsnd * KernelSph::ksrh / this->ksr);
    }

    void calcAlphaDot() {
        PS::F64 src;
        src = - divv * (alphamax - this->alph);
        src = (src > 0.) ? src : 0.;
        this->adot = - (this->alph - alphamin)
            * (0.25 * this->vsnd * KernelSph::ksrh) / this->ksr + src;
    }

    PS::F64 calcTimeStep() {
        return tceff * 2. * this->ksr / (this->vsmx * KernelSph::ksrh);
    }

//    PS::F64 calcTimeStep() {
//        PS::F64 dthydro  = tceff * 2. * this->ksr / this->vsmx;
//        PS::F64 dtenergy = tceff * this->uene / fabs(this->udot);
//        return ((dthydro < dtenergy) ? dthydro : dtenergy);
//    }

    PS::F64 calcEnergy() {
        return this->mass * (0.5 * this->vel * this->vel + this->uene + 0.5 * this->pot);
    }

    static inline PS::F64 calcVolumeInverse(const PS::F64 hi);
    static inline PS::F64 calcPowerOfDimInverse(PS::F64 mass,
                                                PS::F64 dens);

#ifdef WD_DAMPING
    void predict(PS::F64 dt) {
        this->pos   = this->pos  +       this->vel  * dt  + 0.5 * this->acc * dt * dt;
        this->vel2  = this->vel  + 0.5 * this->acc  * dt;
        this->vel   = this->vel  +       this->acc  * dt;
        this->uene2 = this->uene + 0.5 * this->udot * dt;
        this->uene  = this->uene +       this->udot * dt;
        this->alph2 = this->alph + 0.5 * this->adot * dt;
        this->alph  = this->alph +       this->adot * dt;
    }

    void correct(PS::F64 dt) {
        this->acc  -= this->vel / (128.d * dt); //// damping force
        this->vel   = this->vel2  + 0.5 * this->acc  * dt;
        this->uene  = this->uene2 + 0.5 * this->udot * dt;
        this->alph  = this->alph2 + 0.5 * this->adot * dt;
        this->uene  = CalcEquationOfState::getEnergy(this->dens, this->uene); //// damping
    }
#else
    void predict(PS::F64 dt) {
        this->pos   = this->pos  +       this->vel  * dt  + 0.5 * this->acc * dt * dt;
        this->vel2  = this->vel  + 0.5 * this->acc  * dt;
        this->vel   = this->vel  +       this->acc  * dt;
        this->uene2 = this->uene + 0.5 * this->udot * dt;
        this->uene  = this->uene +       this->udot * dt;
        this->alph2 = this->alph + 0.5 * this->adot * dt;
        this->alph  = this->alph +       this->adot * dt;
    }

    void correct(PS::F64 dt) {
        this->vel   = this->vel2  + 0.5 * this->acc  * dt;
        this->uene  = this->uene2 + 0.5 * this->udot * dt;
        this->alph  = this->alph2 + 0.5 * this->adot * dt;
    }
#endif

    inline void dampVelocity(PS::F64 dt) {
        this->vel *= exp(- 0.1 * this->vsnd * KernelSph::ksrh / this->ksr * dt);
    }

};

PS::F64ort SPH::cbox;
PS::F64    SPH::cinv;
PS::F64    SPH::tceff;
PS::F64    SPH::alphamax, SPH::alphamin;
PS::F64    SPH::eps;

#ifdef USE_AT1D
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return mass / dens;
}
#else
#ifdef USE_AT2D
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi * hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return sqrt(mass / dens);
}
#else
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi * hi * hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return pow(mass / dens, 1.d / 3.d);
}
#endif
#endif

template <class Theader>
void setParameterParticle(Theader & header) {
    SPH::cbox     = header.cbox;
    SPH::alphamax = header.alphamax;
    SPH::alphamin = header.alphamin;
    SPH::tceff    = header.tceff;
    SPH::eps      = header.eps;
    
    return;
}
