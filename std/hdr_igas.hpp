#pragma once

namespace CodeUnit {
    PS::F64 UnitOfLength    = 1.;
    PS::F64 UnitOfMass      = 1.;
    PS::F64 UnitOfTime      = 1.;
    PS::F64 GravityConstant = 1.;
}

class CalcEquationOfState {
private:
    CalcEquationOfState() {};
    ~CalcEquationOfState() {};
    CalcEquationOfState(const CalcEquationOfState & c);
    CalcEquationOfState & operator=(const CalcEquationOfState & c);
    static CalcEquationOfState & getInstance() {
        static CalcEquationOfState inst;
        return inst;
    }
public:
    static PS::F64 getEnergyMin(PS::F64 density,
                                PS::F64 abar,
                                PS::F64 zbar) {
        return 0.;
    }
};

class IdealGas : public SPH {
public:
    void readAscii(FILE * fp) {
        fscanf(fp, "%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->mass,
               &this->pos[0], &this->pos[1], &this->pos[2],
               &this->vel[0], &this->vel[1], &this->vel[2],
               &this->uene,   &this->ksr);        
        this->alph  = RP::AlphaInit;
        this->alphu = RP::AlphuInit;
    }

    void writeAscii(FILE *fp) const {
        fprintf(fp, "%6d %+e", this->id, this->mass); // 2
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]); // 5
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]); // 8
        fprintf(fp, " %+e %+e %+e", this->acc[0], this->acc[1], this->acc[2]); // 11
        fprintf(fp, " %+e %+e %+e", this->uene,  this->alph, this->alphu); // 14
        fprintf(fp, " %+e %+e %+e", this->dens, this->vsnd, this->pres); // 17
        fprintf(fp, " %+e %+e %+e", this->divv, this->rotv, this->bswt); // 20
        fprintf(fp, " %+e %+e %+e", this->grdh, this->pot, this->ksr); // 23
        fprintf(fp, " %4d", this->np); // 24
        fprintf(fp, " %+e", this->diffu);
        //fprintf(fp, " %+e %+e", this->vsmx, this->udot);
        fprintf(fp, "\n");
    }

    void referEquationOfState() {
        this->pres  = (RP::AdiabaticIndex - 1.) * this->dens * this->uene;
        this->vsnd  = sqrt(RP::AdiabaticIndex * this->pres / this->dens);
    }

    PS::F64 calcTimeStep() {
        PS::F64 tceff   = RP::CoefficientOfTimestep;
        PS::F64 dth = tceff * this->ksr / (this->vsmx * SK::ksrh);
        PS::F64 dtu = tceff * fabs(this->uene / this->udot);
        return std::min(dth, dtu);
    }

    PS::F64 calcEnergy() {
        return this->mass * (0.5 * this->vel * this->vel + this->uene + 0.5 * this->pot);
    }

    PS::F64vec calcMomentum() {
        return this->mass * this->vel;
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
        FILE *fp = fopen(filename, "r");
        fscanf(fp, "%lf%lf%lf", &Time, &TimeEnd, &TimestepAscii);
        fscanf(fp, "%lld%lld%lld", &NumberOfDimension, &KernelType, &FlagGravity);
        fscanf(fp, "%lf", &ComputationalBox.low_[0]);
        fscanf(fp, "%lf", &ComputationalBox.low_[1]);
        fscanf(fp, "%lf", &ComputationalBox.low_[2]);
        fscanf(fp, "%lf", &ComputationalBox.high_[0]);
        fscanf(fp, "%lf", &ComputationalBox.high_[1]);
        fscanf(fp, "%lf", &ComputationalBox.high_[2]);
        fscanf(fp, "%lf", &AdiabaticIndex);
        fscanf(fp, "%lf%lf%lf", &AlphaMaximum, &AlphaMinimum, &AlphaInit);
        fscanf(fp, "%lf%lf%lf", &AlphuMaximum, &AlphuMinimum, &AlphuInit);
        fscanf(fp, "%lf", &CoefficientOfTimestep);
        fscanf(fp, "%lld", &NumberOfParticle);
        fclose(fp);
    }

    void outputRunParameter(char **argv) {
        FILE * fp = FilePointerForLog;
        fprintf(fp, "# # of process %8d\n", PS::Comm::getNumberOfProc());
        fprintf(fp, "# # of thread  %8d\n", PS::Comm::getNumberOfThread());
        fprintf(fp, "# Header file: %s\n", argv[2]);
        fprintf(fp, "# Data   file: %s\n", argv[3]);
        fprintf(fp, "# KernelType: %d\n", KernelType);
        fprintf(fp, "# AlphaMax:   %+e\n", AlphaMaximum);
        fprintf(fp, "# AlphaMin:   %+e\n", AlphaMinimum);
        fprintf(fp, "# AlphuMax:   %+e\n", AlphuMaximum);
        fprintf(fp, "# AlphuMin:   %+e\n", AlphuMinimum);
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

template <class Tsph>
void outputData(Tsph & sph) {
    if(RP::Time >= RP::TimeAscii) {
        char filename[64];
        sprintf(filename, "snap/sph_t%04d.dat", RP::NumberOfAscii);
        sph.writeParticleAscii(filename);
        RP::TimeAscii += RP::TimestepAscii;
        RP::NumberOfAscii++;        
    }
    PS::F64    etot = calcEnergy(sph);
    PS::F64vec vtot = calcMomentum(sph);
    if(PS::Comm::getRank() == 0) {
        fprintf(RP::FilePointerForLog, "time: %.10f %+e %+e %+e %+e\n",
                RP::Time, RP::Timestep, etot, vtot[0], WT::getTimeTotal());
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
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].correct(RP::Timestep);
    }
}

void initializeSimulation() {
    RP::MaximumTimestep = 1. / 64.;
    RP::MinimumTimestep = 1e-16;
    RP::TimeAscii       = 0.;
    RP::NumberOfStep    = 0;
    RP::NumberOfAscii   = 0;
    RP::NumberOfHexa    = 0;
    //RP::KernelSupportRadiusMaximum;
    //RP::NuclearEnergyTotal;
    RP::Timestep        = RP::MaximumTimestep;
    RP::FilePointerForLog  = fopen("snap/time.log", "w");
    RP::FilePointerForTime = fopen("snap/prof.log", "w");
}

template <class Tdinfo,
          class Tsph,
          class Tdensity,
          class Thydro,
          class Tgravity>
void startSimulation(char **argv,
                     Tdinfo & dinfo,
                     Tsph & sph,
                     Tdensity & density,
                     Thydro & hydro,
                     Tgravity & gravity) {
    RP::readAscii(argv[2]);
    RP::outputRunParameter(argv);
    RP::KernelSupportRadiusMaximum = std::numeric_limits<PS::F64>::max();
    RP::EpsilonOfInternalEnergy    = RP::setEpsilonOfInternalEnergy(sph);
    sph.readParticleAscii(argv[3]);    
    if(RP::ComputationalBox.low_[0] != RP::ComputationalBox.high_[0]) {
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
        dinfo.setPosRootDomain(RP::ComputationalBox.low_, RP::ComputationalBox.high_);
    }
    ND::setDimension(RP::NumberOfDimension);
    SK::setKernel(RP::KernelType, RP::NumberOfDimension);
    sph.adjustPositionIntoRootDomain(dinfo);
    PS::MT::init_genrand(0);
    dinfo.decomposeDomainAll(sph);
    sph.exchangeParticle(dinfo);
    calcSPHKernel(dinfo, sph, density, hydro, gravity);
}

template <class Tdinfo,
          class Tsph>
void restartSimulation(char **argv,
                       Tdinfo & dinfo,
                       Tsph & sph) {
    return;
}

template <class Tdinfo,
          class Tsph,
          class Tdensity,
          class Thydro,
          class Tgravity>
void loopSimulation(Tdinfo & dinfo,
                    Tsph & sph,
                    Tdensity & density,
                    Thydro & hydro,
                    Tgravity & gravity) {

    WT::clear();
    while(RP::Time < RP::TimeEnd) {
        RP::NumberOfStep++;
        RP::Timestep = calcTimeStep(sph);
        WT::reduceInterProcess();
        outputData(sph);
        //PS::Finalize();
        //exit(0);
        WT::clear();
        WT::start();
        predict(sph);
        if(RP::ComputationalBox.low_[0] != RP::ComputationalBox.high_[0]) {
            sph.adjustPositionIntoRootDomain(dinfo);
        }
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
        WT::start();
        calcSPHKernel(dinfo, sph, density, hydro, gravity);
        WT::accumulateCalcHydro();
        WT::start();
        correct(sph);
        WT::accumulateOthers();
        RP::Time += RP::Timestep;
    }
}

template <class Tsph>
void finalizeSimulation(Tsph & sph) {
    outputData(sph);
    fclose(RP::FilePointerForLog);
    fclose(RP::FilePointerForTime);
}

typedef IdealGas GeneralSPH;
