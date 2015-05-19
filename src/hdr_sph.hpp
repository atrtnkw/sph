#pragma once

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
    PS::F64    gamma;    
    PS::F64    alphamax;
    PS::F64    alphamin;
    PS::F64    eta;
    PS::F64    tceff;
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
        eta        = 0.0;
        tceff      = 0.0;
        nptcl      = 0;
    }

    PS::S32 readAscii(FILE *fp) {
        fscanf(fp, "%lf%lf%lf", &time, &tend, &dtsp);
        fscanf(fp, "%lf%lf%lf", &cbox.low_[0], &cbox.low_[1], &cbox.low_[2]);
        fscanf(fp, "%lf%lf%lf", &cbox.high_[0], &cbox.high_[1], &cbox.high_[2]);
        fscanf(fp, "%lf%lf%lf%lf%lf", &gamma, &alphamax, &alphamin, &eta, &tceff);
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
    static PS::F64ort cbox;
    static PS::F64    cinv;
    static PS::F64    ksrh;
    static PS::F64    eta;
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
        this->acc  += gravity.acc;
        this->pot   = gravity.pot + this->mass / this->eps;
        this->accg  = gravity.acc;
    }

    void readAscii(FILE *fp) {
        fscanf(fp, "%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->mass,
               &this->pos[0], &this->pos[1], &this->pos[2],
               &this->vel[0], &this->vel[1], &this->vel[2],
               &this->uene,   &this->alph,   &this->ksr);
    }

    void writeAscii(FILE *fp) const {
        fprintf(fp, "%6d %+e", this->id, this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->acc[0], this->acc[1], this->acc[2]);
        fprintf(fp, " %+e %+e %+e", this->uene, this->alph, this->ksr);
        fprintf(fp, " %+e %+e %+e", this->dens, this->vsnd, this->pres);
        fprintf(fp, " %+e %+e %+e", this->divv, this->rotv, this->bswt);
        fprintf(fp, " %+e %6d %+e", this->grdh, this->np, this->pot);
        fprintf(fp, " %+e %+e %+e", this->accg[0], this->accg[1], this->accg[2]);
        fprintf(fp, "\n");
    }

    void referEquationOfState() {
        this->pres = cinv * this->dens * this->uene;
        this->vsnd = sqrt((1. + cinv) * this->pres / this->dens);
    }

    void calcBalsaraSwitch() {
        this->bswt = fabs(this->grdh * this->divv) / (fabs(this->grdh * this->divv) + fabs(this->grdh * this->rotv) + 1e-4 * this->vsnd / this->ksr);
    }

    void calcAlphaDot() {
        PS::F64 src;
        src = - divv * (alphamax - this->alph);
        src = (src > 0.) ? src : 0.;
        this->adot = - (this->alph - alphamin) * (0.1 * this->vsnd) / this->ksr + src;
    }

    PS::F64 calcTimeStep() {
        return tceff * 2. * this->ksr / this->vsmx;
    }
    /*
    PS::F64 calcTimeStep() {
        PS::F64 dthydro  = tceff * 2. * this->ksr / this->vsmx;
        PS::F64 dtenergy = tceff * this->uene / fabs(this->udot);
        return ((dthydro < dtenergy) ? dthydro : dtenergy);
    }
    */

    PS::F64 calcEnergy() {
#ifdef GRAVITY
        /*
        return this->mass * (0.5 * this->vel * this->vel + this->uene
                             + 0.5 * (this->pot + this->mass / this->eps));
        */
        return this->mass * (0.5 * this->vel * this->vel + this->uene + 0.5 * this->pot);
#else
        return this->mass * (0.5 * this->vel * this->vel + this->uene);
#endif
    }

    static inline PS::F64 calcPowerOfDimInverse(PS::F64 mass,
                                                PS::F64 dens);

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

};

PS::F64ort SPH::cbox;
PS::F64    SPH::cinv;
PS::F64    SPH::eta;
PS::F64    SPH::tceff;
PS::F64    SPH::alphamax, SPH::alphamin;
PS::F64    SPH::eps = 1e-3;

#ifdef USE_AT1D
PS::F64    SPH::ksrh = 1.620185d;
PS::F64    SPH::calcPowerOfDimInverse(PS::F64 mass,
                                      PS::F64 dens) {
    return mass / dens;
}
#else
#ifdef USE_AT2D
PS::F64    SPH::ksrh = 1.897367d;
PS::F64    SPH::calcPowerOfDimInverse(PS::F64 mass,
                                      PS::F64 dens) {
    return sqrt(mass / dens);
}
#else
PS::F64    SPH::ksrh = 1.936492d;
PS::F64    SPH::calcPowerOfDimInverse(PS::F64 mass,
                                      PS::F64 dens) {
    return pow(mass / dens, 1.d / 3.d);
}
#endif
#endif
