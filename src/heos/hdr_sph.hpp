#pragma once

#include "hdr_eos.hpp"

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
    PS::S64    istar;
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
    PS::S32    cnteos;
    static PS::F64    abar;
    static PS::F64    zbar;
    static PS::F64ort cbox;
    static PS::F64    cinv;
    static PS::F64    alphamax, alphamin;
    static PS::F64    tceff;
    static PS::F64    eps;
    static PS::F64vec omg;

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

        fscanf(fp, "%lld%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->istar, &this->mass,
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

        fprintf(fp, "%6d %2d %+e", this->id, this->istar, tmass);
        fprintf(fp, " %+.16e %+.16e %+.16e", tpos[0], tpos[1], tpos[2]);
        fprintf(fp, " %+.16e %+.16e %+.16e", tvel[0], tvel[1], tvel[2]);
        fprintf(fp, " %+.16e %+.16e %+.16e", tacc[0], tacc[1], tacc[2]);
        fprintf(fp, " %+.16e %+.16e %+.16e", tuene, this->alph, tksr);
        fprintf(fp, " %+.16e %+.16e %+.16e %+.16e", tdens, tvsnd, tpres, this->temp);
        fprintf(fp, " %+.16e %+.16e %+.16e", tdivv, trotv, this->bswt);
        fprintf(fp, " %+.16e %6d %+.16e", this->grdh, this->np, tpot);        
//        PS::F64vec taccg = this->accg * UnitOfAcceleration;
//        fprintf(fp, " %+.16e %+.16e %+.16e", taccg[0], taccg[1], taccg[2]);
        fprintf(fp, " %3d", this->cnteos);
        fprintf(fp, "\n");

    }

    void referEquationOfState() {
#ifdef WD_DAMPING1
        this->pres = CalcEquationOfState::getPressure(this->dens, this->uene);
        this->vsnd = CalcEquationOfState::getSoundVelocity(this->dens, this->uene);
        this->temp = CalcEquationOfState::getTemperature(this->dens, this->uene);
#else
        CalcEquationOfState::getThermodynamicQuantity(this->dens,
                                                      this->uene,
                                                      this->abar,
                                                      this->zbar,
                                                      this->pres,
                                                      this->vsnd,
                                                      this->temp,
                                                      this->cnteos);
#endif
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

#ifdef WD_DAMPINGB
    inline void addAdditionalForce() {
        this->acc  -= this->omg ^ (this->omg ^ this->pos) + 2.d * (this->omg ^ this->vel);
    }

    PS::F64 calcEnergy() {
        PS::F64vec tv = this->omg ^ this->pos;
        return this->mass * (0.5 * this->vel * this->vel + this->uene
                             + 0.5 * (this->pot - tv * tv));
    }
#else
    inline void addAdditionalForce() {
       ;
    }
    
    PS::F64 calcEnergy() {
        return this->mass * (0.5 * this->vel * this->vel + this->uene + 0.5 * this->pot);
    }
#endif

    static inline PS::F64 calcVolumeInverse(const PS::F64 hi);
    static inline PS::F64 calcPowerOfDimInverse(PS::F64 mass,
                                                PS::F64 dens);
    static inline v4df calcVolumeInverse(const v4df hi);

    void predict(PS::F64 dt) {
        this->pos   = this->pos  +       this->vel  * dt  + 0.5 * this->acc * dt * dt;
        this->vel2  = this->vel  + 0.5 * this->acc  * dt;
        this->vel   = this->vel  +       this->acc  * dt;
        this->uene2 = this->uene + 0.5 * this->udot * dt;
        this->uene  = this->uene +       this->udot * dt;
        this->alph2 = this->alph + 0.5 * this->adot * dt;
        this->alph  = this->alph +       this->adot * dt;
    }
    
#ifdef WD_DAMPING1
    void correct(PS::F64 dt) {
        this->acc  -= this->vel / (128.d * dt);
        this->vel   = this->vel2  + 0.5 * this->acc  * dt;
        this->uene  = this->uene2 + 0.5 * this->udot * dt;
        this->alph  = this->alph2 + 0.5 * this->adot * dt;
        this->uene  = CalcEquationOfState::getEnergy(this->dens, this->uene);
    }
#elif defined WD_DAMPING2
    void correct(PS::F64 dt) {
        this->vel   = this->vel2  + 0.5 * this->acc  * dt;
        this->uene  = this->uene2 + 0.5 * this->udot * dt;
        this->alph  = this->alph2 + 0.5 * this->adot * dt;

        PS::F64 unow = this->uene;
        PS::F64 umin = CalcEquationOfState::getEnergyMin(this->dens, this->abar, this->zbar);
        PS::F64 delu = (unow < umin) ? 0.d : unow - umin;
        this->uene = (unow < umin) ? unow : ((unow - umin) * exp(-0.1 * dt) + umin);
    }
#elif defined WD_DAMPINGB
    void correct(PS::F64 dt) {
        this->acc  -= this->vel / (128.d * dt);
        this->vel   = this->vel2  + 0.5 * this->acc  * dt;
        this->uene  = this->uene2 + 0.5 * this->udot * dt;
        this->alph  = this->alph2 + 0.5 * this->adot * dt;
    }
#else
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

PS::F64    SPH::abar = 13.7142857143d;
PS::F64    SPH::zbar =  6.85714285714d;
PS::F64ort SPH::cbox;
PS::F64    SPH::cinv;
PS::F64    SPH::tceff;
PS::F64    SPH::alphamax, SPH::alphamin;
PS::F64    SPH::eps;
PS::F64vec SPH::omg;

#ifdef USE_AT1D
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return mass / dens;
}
inline v4df SPH::calcVolumeInverse(const v4df hi) {return hi;}
#else
#ifdef USE_AT2D
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi * hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return sqrt(mass / dens);
}
inline v4df SPH::calcVolumeInverse(const v4df hi) {return hi * hi;}
#else
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi * hi * hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return pow(mass / dens, 1.d / 3.d);
}
inline v4df SPH::calcVolumeInverse(const v4df hi) {return hi * hi * hi;}
#endif
#endif

template <class Tptcl>
void calcCenterOfMass(Tptcl & system,
                      PS::F64    & mc,
                      PS::F64vec & xc,
                      PS::F64vec & vc,
                      PS::S64 istar = -1) {
    PS::F64    mloc = 0.d;
    PS::F64vec xloc = 0.d;
    PS::F64vec vloc = 0.d;

    assert(istar < 2);

    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        if(system[i].istar == istar || istar < 0) {
            mloc += system[i].mass;
            xloc += system[i].mass * system[i].pos;
            vloc += system[i].mass * system[i].vel;
        }
    }

    mc = PS::Comm::getSum(mloc);
    xc = PS::Comm::getSum(xloc);
    vc = PS::Comm::getSum(vloc);
    PS::F64 minv = 1.d / mc;
    xc *= minv;
    vc *= minv;    
}

template <class Theader>
void setParameterParticle(Theader & header) {
    SPH::cbox     = header.cbox;
    SPH::alphamax = header.alphamax;
    SPH::alphamin = header.alphamin;
    SPH::tceff    = header.tceff;
    SPH::eps      = header.eps;
    
    return;
}

template <class Tptcl>
void reduceSeparation(PS::F64 time, 
                      Tptcl & system,
                      FILE * fplog = NULL) {

#ifdef WD_DAMPINGB
    static bool    firststep = true;
    static PS::F64 ReductionTime;
    static PS::F64 DeltaSystemTime = 1.d / 64.d;
//    static PS::F64 CriticalRadius  = 1.5e9 * CodeUnit::UnitOfLengthInv;
    static PS::F64 CriticalRadius  = 1.8e9 * CodeUnit::UnitOfLengthInv;
    static bool StopDampingB = false;

    if(firststep) {
        PS::F64    m1;
        PS::F64vec x1;
        PS::F64vec v1;
        calcCenterOfMass(system, m1, x1, v1, 1);

        PS::F64 r2loc = 0.d;
        PS::S32 nloc = system.getNumberOfParticleLocal();
        for(PS::S32 i = 0; i < nloc; i++) {
            if(system[i].istar == 1) {
                PS::F64vec dr = system[i].pos - x1;
                PS::F64    r2 = dr * dr;
                if(r2 > r2loc) {
                    r2loc = r2;
                }
            }
        }        

        PS::F64 r2glb = PS::Comm::getMaxValue(r2loc);
        PS::F64 rmax  = sqrt(r2glb);
        PS::F64 rho2  = m1 / ((4.d * M_PI / 3.d * rmax * rmax * rmax));

        ReductionTime = 1.d / (0.05 * sqrt(CodeUnit::grav * rho2));;

        firststep = false;
    }

    if(time / DeltaSystemTime - (PS::S64)(time / DeltaSystemTime) == 0.d
        && time != 0.d) {
        PS::F64    m0, m1;
        PS::F64vec x0, x1;
        PS::F64vec v0, v1;
        calcCenterOfMass(system, m0, x0, v0, 0);
        calcCenterOfMass(system, m1, x1, v1, 1);

        if(StopDampingB) {
            PS::F64 mc;
            PS::F64vec xc;
            PS::F64vec vc;
            calcCenterOfMass(system, mc, xc, vc);            
            PS::S32    nloc  = system.getNumberOfParticleLocal();
            for(PS::S32 i = 0; i < nloc; i++) {
                system[i].pos -= xc;
                system[i].vel -= vc;
            }    

            calcCenterOfMass(system, m0, x0, v0, 0);
            calcCenterOfMass(system, m1, x1, v1, 1);
            
            for(PS::S32 i = 0; i < nloc; i++) {                
                system[i].vel += system[i].omg ^ system[i].pos;
            }

            char filename[64];
            sprintf(filename, "snap/final.dat");
            system.writeParticleAscii(filename);

            PS::Finalize();
            exit(0);
        }
    
        PS::F64vec dx = x1 - x0;
        PS::F64    dr = sqrt(dx * dx);
        PS::F64    ds = dr / ReductionTime * DeltaSystemTime;
    
        dx *= (- ds / dr);

        PS::S32 nloc = system.getNumberOfParticleLocal();
        for(PS::S32 i = 0; i < nloc; i++) {
            if(system[i].istar == 1) {
                system[i].pos += dx;
            }
        }

        PS::F64    mc;
        PS::F64vec xc;
        PS::F64vec vc;
        calcCenterOfMass(system, mc, xc, vc);

        for(PS::S32 i = 0; i < nloc; i++) {
            system[i].pos -= xc;
        }

        calcCenterOfMass(system, m0, x0, v0, 0);
        calcCenterOfMass(system, m1, x1, v1, 1);
        dx = x1 - x0;
        dr = sqrt(dx * dx);

        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "bsep: %.10f %+e\n", time, dr * CodeUnit::UnitOfLength);
            fprintf(fplog,  "bsep: %.10f %+e\n", time, dr * CodeUnit::UnitOfLength);
            fflush(fplog);
        }

        if(dr < CriticalRadius) {
            StopDampingB = true;
        }
   }
#endif    
    
}

template <class Tptcl>
void calcFieldVariable(Tptcl & system) {

#ifdef WD_DAMPINGB
    PS::F64    m0, m1;
    PS::F64vec x0, x1;
    PS::F64vec v0, v1;
    calcCenterOfMass(system, m0, x0, v0, 0);
    calcCenterOfMass(system, m1, x1, v1, 1);

    PS::F64vec axisv = x0 - x1;
    PS::F64    axis  = sqrt(axisv * axisv);
    PS::F64    vel   = sqrt(CodeUnit::grav * (m0 + m1) / axis);
    
    SPH::omg[0] = 0.d;
    SPH::omg[1] = 0.d;
    SPH::omg[2] = vel / axis;

#endif

    return;
}

template <class Tptcl>
void doThisEveryTime(PS::F64 & time,
                     PS::F64 & dtime,
                     PS::F64 & tout,
                     PS::F64 & dtsp,
                     Tptcl & system,
                     FILE * fplog,
                     FILE * fptim) {
    reduceSeparation(time, system, fplog);

    calcFieldVariable(system);

    if(time >= tout) {
        char filename[64];
        sprintf(filename, "snap/sph_t%04d.dat", (PS::S32)time);
        system.writeParticleAscii(filename);
        tout += dtsp;
    }

    if(time == 1.0) {
        PS::Finalize();
        exit(0);
    }
    
    PS::F64 etot = calcEnergy(system);
    WT::reduceInterProcess();
    if(PS::Comm::getRank() == 0) {
        using namespace CodeUnit;
        fprintf(fplog,  "time: %.10f %+e %+e\n", time * UnitOfTime,
                etot * UnitOfEnergy * UnitOfMass, WT::getTimeTotal());
        fflush(fplog);
        WT::dump(time, fptim);
        fflush(fptim);
    }
    WT::clear();

}

template <class Tptcl>
void finalizeSimulation(PS::S32 time,
                        Tptcl & system) {
    char filename[64];

    sprintf(filename, "snap/sph_t%04d.dat", (PS::S32)time);
    system.writeParticleAscii(filename);

    PS::F64    mc;
    PS::F64vec xc;
    PS::F64vec vc;
    calcCenterOfMass(system, mc, xc, vc);
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        system[i].pos -= xc;
        system[i].vel -= vc;
    }

    sprintf(filename, "snap/final.dat");
    system.writeParticleAscii(filename);
}

