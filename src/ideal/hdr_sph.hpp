#pragma once

namespace CodeUnit {
    PS::F64 grav         = 1.;
    PS::F64 UnitOfLength = 1.;
    PS::F64 UnitOfEnergy = 1.;
}

class CalcEquationOfState {
public:
    static PS::F64 getEnergyMin(PS::F64 dens,
                         PS::F64 abar,
                         PS::F64 zbar) {
        return 0.;
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

class Gradient {
public:
    PS::F64mat ctau;
    void clear() {
        ctau = 0.;
    }
};

class Derivative{
public:
    PS::F64vec acc;
    PS::F64vec accg;
    PS::F64    udot;
    PS::F64    vsmx;
    PS::F64    diffu;
    void clear(){
        acc   = 0.0;
        accg  = 0.0;
        udot  = 0.0;
        vsmx  = 0.0;
        diffu = 0.0;
    }
};

class Gravity{
public:
    PS::F64vec acc;
    PS::F64    pot;
    PS::F64    eta;
    void clear(){
        acc = 0.0;
        pot = 0.0;
        eta = 0.0;
    }
};

class Header {
public:
    PS::F64    time;
    PS::F64    tend;
    PS::F64    dtsp;
    PS::F64ort cbox;
    PS::F64    gamma;    
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
        gamma      = 0.0;
        alphamax   = 0.0;
        alphamin   = 0.0;
        tceff      = 0.0;
        eps        = 0.0;
        nptcl      = 0;
    }

    PS::S32 readAscii(FILE *fp) {
        fscanf(fp, "%lf%lf%lf", &time, &tend, &dtsp);
        fscanf(fp, "%lf%lf%lf", &cbox.low_[0], &cbox.low_[1], &cbox.low_[2]);
        fscanf(fp, "%lf%lf%lf", &cbox.high_[0], &cbox.high_[1], &cbox.high_[2]);
        fscanf(fp, "%lf%lf%lf%lf", &gamma, &alphamax, &alphamin, &tceff);
#ifdef GRAVITY
        fscanf(fp, "%lf", &eps);
#endif
        fscanf(fp, "%d", &nptcl);

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
    PS::F64    alphu;
    PS::F64    alphu2;
    PS::F64    adotu;
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
    PS::F64vec acch;
    PS::F64vec accg1;
    PS::F64vec accg2;
    PS::F64    eta;
    PS::F64    diffu;
    PS::F64    abar, zbar;
    PS::F64mat ctau;
    static PS::F64ort cbox;
    static PS::F64    cinv;
    static PS::F64    alphamax, alphamin;
    static PS::F64    tceff;
    static PS::F64    eps;
    static PS::F64    ksrmax;
    static PS::F64vec omg;
    static PS::F64    epsu;

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

    void copyFromForce(const Gradient & gradient) {
        this->ctau = gradient.ctau;
    }

    void copyFromForce(const Gravity & gravity) {
        this->accg2 = gravity.acc;
        this->eta   = this->ksr * this->ksr * KernelSph::ksrhinv * KernelSph::ksrhinv
            * this->grdh / (KernelSph::dim * this->dens) * gravity.eta;
        this->pot   = gravity.pot;
    }

    void copyFromForce(const Derivative & derivative){
        this->acch  = derivative.acc;
        this->accg1 = derivative.accg;
        this->udot  = derivative.udot;
        this->vsmx  = derivative.vsmx;
        this->acc   = this->acch + this->accg1 + this->accg2;
        this->diffu = derivative.diffu;
    }

    void readAscii(FILE *fp) {
        fscanf(fp, "%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->mass,
               &this->pos[0], &this->pos[1], &this->pos[2],
               &this->vel[0], &this->vel[1], &this->vel[2],
               &this->uene,   &this->alph,   &this->ksr);
        this->alphu = this->alph;
    }

    void writeAscii(FILE *fp) const {
        fprintf(fp, "%6d %+e", this->id, this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->acc[0], this->acc[1], this->acc[2]);
        fprintf(fp, " %+e %+e %+e %+e", this->uene, this->alph, this->alphu, this->ksr);
        fprintf(fp, " %+e %+e %+e", this->dens, this->vsnd, this->pres);
        fprintf(fp, " %+e %+e %+e", this->divv, this->rotv, this->bswt);
        fprintf(fp, " %+e %6d %+e", this->grdh, this->np, this->pot);
        fprintf(fp, " %+e", this->eta);
        /////////////////////
        fprintf(fp, " %+e %+e %+e", this->adotu, this->diffu, this->udot);
        fprintf(fp, " %+e %+e %+e", this->ctau.xx, this->ctau.yy, this->ctau.zz);
        fprintf(fp, " %+e %+e %+e", this->ctau.xy, this->ctau.xz, this->ctau.yz);
        /////////////////////
        fprintf(fp, "\n");
    }

    void referEquationOfState() {
        this->pres = cinv * this->dens * this->uene;
        this->vsnd = sqrt((1. + cinv) * this->pres / this->dens);
    }

    void calcBalsaraSwitch() {
        this->bswt = fabs(this->divv)
            / (fabs(this->divv) + fabs(this->rotv)
               + 1e-4 * this->vsnd * KernelSph::ksrh / this->ksr);
    }

    void calcAlphaDot() {
        PS::F64 src;
        src = - divv * (alphamax - this->alph);
        src = (src > 0.) ? src : 0.;
        PS::F64 tauinv = (0.25 * this->vsnd * KernelSph::ksrh) / this->ksr;
        this->adot = - (this->alph - alphamin) * tauinv + src;
#ifdef THERMAL_CONDUCTIVITY
        PS::F64 srcu;
        srcu  = this->ksr * KernelSph::ksrhinv * this->diffu / sqrt(this->uene + this->epsu);
        srcu *= (alphamax - this->alphu);
        this->adotu = - (this->alphu - alphamin) * tauinv + srcu;
#endif
    }

    // check
    PS::F64 calcTimeStep() {
        PS::F64 dthydro = tceff * 2. * this->ksr / (this->vsmx * KernelSph::ksrh);
        PS::F64 dtenerg = tceff * this->uene / fabs(this->udot);
        return std::min(dthydro, dtenerg);
    }

    PS::F64 calcEnergy() {
        return this->mass * (0.5 * this->vel * this->vel + this->uene + 0.5 * this->pot);
    }

    static inline PS::F64 calcVolumeInverse(const PS::F64 hi);
    static inline PS::F64 calcPowerOfDimInverse(PS::F64 mass,
                                                PS::F64 dens);
    static inline v4df calcVolumeInverse(const v4df hi);
    static inline v8sf calcVolumeInverse(const v8sf hi);

#ifdef DAMPING
    inline void addAdditionalForce() {
        ;
    }

    void predict(PS::F64 dt) {
        this->pos   = this->pos  +       this->vel  * dt  + 0.5 * this->acc * dt * dt;
        this->vel2  = this->vel  + 0.5 * this->acc  * dt;
        this->vel   = this->vel  +       this->acc  * dt;
    }

    void correct(PS::F64 dt) {
        this->vel   = this->vel2  + 0.5 * this->acc  * dt;
    }
#else
    inline void addAdditionalForce() {
        ;
    }

    void predict(PS::F64 dt) {
        this->pos    = this->pos   +       this->vel   * dt  + 0.5 * this->acc * dt * dt;
        this->vel2   = this->vel   + 0.5 * this->acc   * dt;
        this->vel    = this->vel   +       this->acc   * dt;
        this->uene2  = this->uene  + 0.5 * this->udot  * dt;
        this->uene   = this->uene  +       this->udot  * dt;
        this->alph2  = this->alph  + 0.5 * this->adot  * dt;
        this->alph   = this->alph  +       this->adot  * dt;
        this->alphu2 = this->alphu + 0.5 * this->adotu * dt;
        this->alphu  = this->alphu +       this->adotu * dt;
    }

    void correct(PS::F64 dt) {
        this->vel   = this->vel2   + 0.5 * this->acc   * dt;
        this->uene  = this->uene2  + 0.5 * this->udot  * dt;
        this->alph  = this->alph2  + 0.5 * this->adot  * dt;
        this->alphu = this->alphu2 + 0.5 * this->adotu * dt;
    }
#endif

    void dampVelocity(PS::F64 dt) {
        this->vel *= exp(- 0.1 * this->vsnd * KernelSph::ksrh / this->ksr * dt);
    }

};

PS::F64ort SPH::cbox;
PS::F64    SPH::cinv;
PS::F64    SPH::tceff;
PS::F64    SPH::alphamax, SPH::alphamin;
PS::F64    SPH::eps;
PS::F64    SPH::ksrmax = std::numeric_limits<double>::max();
PS::F64vec SPH::omg;
PS::F64    SPH::epsu;

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

template <class Theader>
void setParameterParticle(Theader & header) {
    SPH::cbox     = header.cbox;
    SPH::cinv     = header.gamma - 1.d;
    SPH::alphamax = header.alphamax;
    SPH::alphamin = header.alphamin;
    SPH::tceff    = header.tceff;
    SPH::eps      = header.eps;
    
    return;
}

template <class Tptcl>
PS::F64 setEpsilonOfInternalEnergy(Tptcl & system) {
    PS::F64 uminloc = std::numeric_limits<double>::max();
    PS::S32 nloc    = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        if(system[i].uene < uminloc) {
            uminloc = system[i].uene;
        }
    }
    PS::F64 uminglb = PS::Comm::getMinValue(uminloc);
    return 1e-4 * uminglb;
}

template <class Tptcl>
void calcFieldVariable(Tptcl & system) {
    return;
}

template <class Tptcl>
void finalizeSimulation(PS::S32 nstp,
                        Tptcl & system) {
    char filename[64];

    sprintf(filename, "snap/sph_t%04d.dat", nstp);
    system.writeParticleAscii(filename);

    PS::F64    msloc = 0.0d;
    PS::F64vec xcloc = 0.0d;    
    PS::F64vec vcloc = 0.0d;    
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        msloc += system[i].mass;
        xcloc += system[i].mass * system[i].pos;
        vcloc += system[i].mass * system[i].vel;
    }
}

