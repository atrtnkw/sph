#pragma once

#include "hdr_eos.hpp"

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
    PS::S32    istp;
    PS::F64ort cbox;
    PS::F64    alphamax;
    PS::F64    alphamin;
    PS::F64    tceff;
    PS::F64    eps;
    PS::S64    nptcl;
    PS::F64    ksrmax;

    Header() {
        time       = 0.0;
        tend       = 0.0;
        dtsp       = 0.0;
        istp       = 0;
        cbox.low_  = 0.0;
        cbox.high_ = 0.0;
        alphamax   = 0.0;
        alphamin   = 0.0;
        tceff      = 0.0;
        eps        = 0.0;
        nptcl      = 0;
        ksrmax     = 0.0;
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

    template <class Tdinfo>
    void readRestartFile(FILE *fp,
                         Tdinfo & dinfo) {
        PS::F64 (*cvt)(PS::U64) = convertU64ToF64;        
        PS::U64 utime, utend, udtsp;
        PS::U64 ualphamax, ualphamin, utceff;
        PS::U64 ueps;
        PS::U64 uksrmax;
        fscanf(fp, "%llx %llx %llx %d", &utime, &utend, &udtsp, &this->istp);
        fscanf(fp, "%llx %llx %llx", &ualphamax, &ualphamin, &utceff);
        fscanf(fp, "%llx", &ueps);
        fscanf(fp, "%llx", &uksrmax);
        fscanf(fp, "%d", &this->nptcl);
        time = cvt(utime);
        tend = cvt(utend);
        dtsp = cvt(udtsp);
        alphamax = cvt(ualphamax);
        alphamin = cvt(ualphamin);
        tceff    = cvt(utceff);
        eps = cvt(ueps);
        ksrmax = cvt(uksrmax);

        PS::S32 nproc = PS::Comm::getNumberOfProc();
        for(PS::S32 i = 0; i < nproc; i++) {
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

    template <class Tdinfo>
    void writeRestartFile(FILE *fp,
                          Tdinfo & dinfo) {
        PS::U64 (*cvt)(PS::F64) = convertF64ToU64;
        fprintf(fp, "%llx %llx %llx %d\n", cvt(time), cvt(tend), cvt(dtsp), istp);
        fprintf(fp, "%llx %llx %llx\n", cvt(alphamax), cvt(alphamin), cvt(tceff));
        fprintf(fp, "%llx\n", cvt(eps));
        fprintf(fp, "%llx\n", cvt(ksrmax));
        fprintf(fp, "%d\n", nptcl);

        PS::S32 nproc = PS::Comm::getNumberOfProc();
        for(PS::S32 i = 0; i < nproc; i++) {
            PS::F64ort domain = dinfo.getPosDomain(i);
            fprintf(fp, "%llx %llx %llx %llx %llx %llx\n",
                    cvt(domain.low_[0]),  cvt(domain.low_[1]),  cvt(domain.low_[2]),
                    cvt(domain.high_[0]), cvt(domain.high_[1]), cvt(domain.high_[2]));
        }
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
    static PS::F64    ReductionTimeInv;
    static PS::F64    ksrmax;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {        
        this->pos = pos_new;
    }

    void copyFromForce(const Density & density){
        this->dens  = density.dens;
        this->grdh  = density.grdh;
        this->np    = density.np;
        this->ksr   = density.ksr;
        this->rotv  = density.rotv;
        this->divv  = density.divv;
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
        fprintf(fp, " %3d", this->cnteos);
#ifdef WD_DAMPINGB
        fprintf(fp, " %+.16e", SPH::omg[2]);
#endif
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
    static inline v8sf calcVolumeInverse(const v8sf hi);

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
        //this->acc  -= this->vel / (128.d * dt);
        this->acc  -= this->vel * SPH::ReductionTimeInv;
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

    void readRestartFile(FILE *fp) {
        PS::F64 (*cvt)(PS::U64) = convertU64ToF64;
        PS::S32 s0, s1;
        PS::U64 umass;
        PS::U64 upos[3], uvel[3], uacc[3];
        PS::U64 uuene, uudot;
        PS::U64 ualph, uadot;
        PS::U64 udens, upres, uvsnd, utemp;
        PS::U64 udivv, urotv, ubswt;
        PS::U64 uksr, ugrdh, uvsmx;
        PS::U64 upot;

        fscanf(fp, "%d %d %llx", &this->id, &this->istar, &umass);
        fscanf(fp, "%llx %llx %llx", &upos[0], &upos[1], &upos[2]);
        fscanf(fp, "%llx %llx %llx", &uvel[0], &uvel[1], &uvel[2]);
        fscanf(fp, "%llx %llx %llx", &uacc[0], &uacc[1], &uacc[2]);
        fscanf(fp, "%llx %llx", &uuene, &uudot);
        fscanf(fp, "%llx %llx", &ualph, &uadot);
        fscanf(fp, "%llx %llx %llx %llx", &udens, &upres, &uvsnd, &utemp);
        fscanf(fp, "%llx %llx %llx", &udivv, &urotv, &ubswt);
        fscanf(fp, "%llx %llx %llx", &uksr,  &ugrdh, &uvsmx);
        fscanf(fp, "%llx", &upot);

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
        this->uene = cvt(uuene);
        this->udot = cvt(uudot);
        this->alph = cvt(ualph);
        this->adot = cvt(uadot);
        this->dens = cvt(udens);
        this->pres = cvt(upres);
        this->vsnd = cvt(uvsnd);
        this->temp = cvt(utemp);
        this->divv = cvt(udivv);
        this->rotv = cvt(urotv);
        this->bswt = cvt(ubswt);
        this->ksr  = cvt(uksr);
        this->grdh = cvt(ugrdh);
        this->vsmx = cvt(uvsmx);
        this->pot  = cvt(upot);        
    }

    void writeRestartFile(FILE *fp) const {
        PS::U64 (*cvt)(PS::F64) = convertF64ToU64;
        fprintf(fp, "%6d %2d %llx", this->id, this->istar, cvt(this->mass));
        fprintf(fp, " %llx %llx %llx", cvt(this->pos[0]), cvt(this->pos[1]), cvt(this->pos[2]));
        fprintf(fp, " %llx %llx %llx", cvt(this->vel[0]), cvt(this->vel[1]), cvt(this->vel[2]));
        fprintf(fp, " %llx %llx %llx", cvt(this->acc[0]), cvt(this->acc[1]), cvt(this->acc[2]));
        fprintf(fp, " %llx %llx", cvt(this->uene), cvt(this->udot));
        fprintf(fp, " %llx %llx", cvt(this->alph), cvt(this->adot));
        fprintf(fp, " %llx %llx", cvt(this->dens), cvt(this->pres));
        fprintf(fp, " %llx %llx", cvt(this->vsnd), cvt(this->temp));
        fprintf(fp, " %llx %llx %llx", cvt(this->divv), cvt(this->rotv), cvt(this->bswt));
        fprintf(fp, " %llx %llx %llx", cvt(this->ksr), cvt(this->grdh), cvt(this->vsmx));
        fprintf(fp, " %llx", cvt(this->pot));
        fprintf(fp, "\n");
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
PS::F64    SPH::ReductionTimeInv;
PS::F64    SPH::ksrmax;

#ifdef USE_AT1D
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return mass / dens;
}
inline v4df SPH::calcVolumeInverse(const v4df hi) {return hi;}
inline v8sf SPH::calcVolumeInverse(const v8sf hi) {return hi;}
#else
#ifdef USE_AT2D
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi * hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return sqrt(mass / dens);
}
inline v4df SPH::calcVolumeInverse(const v4df hi) {return hi * hi;}
inline v8sf SPH::calcVolumeInverse(const v8sf hi) {return hi * hi;}
#else
inline PS::F64 SPH::calcVolumeInverse(const PS::F64 hi) {return hi * hi * hi;}
inline PS::F64 SPH::calcPowerOfDimInverse(PS::F64 mass,
                                          PS::F64 dens) {
    return pow(mass / dens, 1.d / 3.d);
}
inline v4df SPH::calcVolumeInverse(const v4df hi) {return hi * hi * hi;}
inline v8sf SPH::calcVolumeInverse(const v8sf hi) {return hi * hi * hi;}
#endif
#endif

template <class Theader,
          class Tdinfo,
          class Tpsys>
void readRestartFile(char * const filename,
                     Theader & header,
                     Tdinfo & dinfo,
                     Tpsys & system) {

    FILE *fp = fopen(filename, "r");

    header.readRestartFile(fp, dinfo);
    PS::S32 nloc = header.nptcl;
    system.setNumberOfParticleLocal(nloc);
    for(PS::S32 i = 0; i < nloc; i++) {
        system[i].readRestartFile(fp);
    }

    fclose(fp);
}

template <class Theader,
          class Tdinfo,
          class Tpsys>
void writeRestartFile(char * const filename,
                      const PS::F64 time,
                      Theader & header,
                      Tdinfo & dinfo,
                      Tpsys & system) {

    FILE *fp = fopen(filename, "w");

    PS::S32 nloc = system.getNumberOfParticleLocal();

    header.time  = time;
    header.nptcl = nloc;

    header.writeRestartFile(fp, dinfo);

    for(PS::S32 i = 0; i < nloc; i++) {
        system[i].writeRestartFile(fp);
    }

    fclose(fp);
}

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
    SPH::ksrmax   = header.ksrmax;
    
    return;
}

template <class Tpsys>
PS::F64 calcSystemSize(Tpsys & system) {
    PS::F64    m0, m1;
    PS::F64vec x0, x1, v0, v1;
    calcCenterOfMass(system, m0, x0, v0, 0);
    calcCenterOfMass(system, m1, x1, v1, 1);
    PS::F64vec dx = x0 - x1;
    PS::F64    dr = (m0 != 0. && m1 != 0.) ? sqrt(dx * dx)
        : std::numeric_limits<double>::max();
    return dr;
}

template <class Tptcl,
          class Tmsls>
void reduceSeparation(PS::F64 time, 
                      Tptcl & system,
                      Tmsls & msls,
                      FILE * fplog = NULL) {

#ifdef WD_DAMPINGB
    static bool    firststep = true;
    static PS::F64 ReductionTime;
    static PS::F64 DeltaSystemTime = MinimumTimeStep;
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
        SPH::ReductionTimeInv = 1. / ReductionTime;

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
            sprintf(filename, "snap/msls_final.dat");
            msls.writeParticleAscii(filename);

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

        PS::F64 xlag = searchLagrange1(msls, x0[0], x1[0]);

        if(PS::Comm::getRank() == 0) {
            fprintf(fplog,  "bsep: %.10f %+e\n", time, dr * CodeUnit::UnitOfLength);
            fflush(fplog);
        }
        
#if 1        
        bool StopDampingBLocal = false;
        for(PS::S32 i = 0; i < nloc; i++) {
            if(system[i].istar == 0)
                continue;
            if(system[i].pos[0] > xlag) {
                StopDampingBLocal = true;
                break;
            }
        }
        StopDampingB = PS::Comm::synchronizeConditionalBranchOR(StopDampingBLocal);
#else
        if(dr < CriticalRadius) {
            StopDampingB = true;
        }
#endif
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

template <class Theader,
          class Tdinfo,
          class Tptcl,
          class Tmassless>
void doThisEveryTime(PS::F64 & dtime,
                     PS::F64 & tout,
                     Theader & header,
                     Tdinfo & dinfo,
                     Tptcl & system,
                     Tmassless & msls,
                     FILE * fplog,
                     FILE * fptim) {
    reduceSeparation(header.time, system, msls, fplog);

    calcFieldVariable(system);

    if(PS::Comm::getRank() == 0 && header.time == 0.) {
        fprintf(fplog, "# Maximum kernel length: %e\n", header.ksrmax);
        fprintf(fplog, "# Unit of length:        %e\n", CodeUnit::UnitOfLength);
    }

    PS::F64 trst = 10.;
    if(header.time - (PS::S64)(header.time / trst) * trst == 0.) {
        char filename[64];
        sprintf(filename, "snap/t%04d_p%06d.hexa", (PS::S32)header.time, PS::Comm::getRank());
        writeRestartFile(filename, header.time, header, dinfo, system);
    }

    if(header.time >= tout) {
        char filename[64];
        sprintf(filename, "snap/sph_t%04d.dat", (PS::S32)header.time);
        system.writeParticleAscii(filename);
#ifdef WD_DAMPINGB
        sprintf(filename, "snap/msls_t%04d.dat", (PS::S32)header.time);
        msls.writeParticleAscii(filename);
#endif
        tout += header.dtsp;
    }

    PS::F64 etot = calcEnergy(system);
    WT::reduceInterProcess();
    if(PS::Comm::getRank() == 0)     {
        using namespace CodeUnit;
        fprintf(fplog,  "time: %.10f %+e %+e\n", header.time * UnitOfTime,
                etot * UnitOfEnergy * UnitOfMass, WT::getTimeTotal());
        fflush(fplog);
        WT::dump(header.time, fptim);
        fflush(fptim);
    }
    WT::clear();

}

template <class Theader,
          class Tdinfo,
          class Tptcl,
          class Tmsls>
void finalizeSimulation(PS::S32 time,
                        Theader & header,
                        Tdinfo & dinfo,
                        Tptcl & system,
                        Tmsls & msls) {

    {
        char filename[64];
        sprintf(filename, "snap/t%04d_p%06d.hexa", (PS::S32)time, PS::Comm::getRank());
        writeRestartFile(filename, time, header, dinfo, system);
    }

    {
        char filename[64];
        sprintf(filename, "snap/sph_t%04d.dat", (PS::S32)time);
        system.writeParticleAscii(filename);
    }

    PS::F64    mc;
    PS::F64vec xc;
    PS::F64vec vc;
    calcCenterOfMass(system, mc, xc, vc);
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        system[i].pos -= xc;
        system[i].vel -= vc;
    }
    {
        char filename[64];
        sprintf(filename, "snap/final.dat");
        system.writeParticleAscii(filename);
    }
}
