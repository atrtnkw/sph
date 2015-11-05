#pragma once

class IdealGas : public SPH {
public:
    void readAscii(FILE * fp) {
        fscanf(fp, "%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->mass,
               &this->pos[0], &this->pos[1], &this->pos[2],
               &this->vel[0], &this->vel[1], &this->vel[2],
               &this->uene,   &this->alph,   &this->ksr);        
    }

    void writeAscii(FILE *fp) const {
        fprintf(fp, "%6d %+e", this->id, this->mass); // 2
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]); // 5
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]); // 8
        fprintf(fp, " %+e %+e %+e", this->acc[0], this->acc[1], this->acc[2]); // 11
        fprintf(fp, " %+e %+e %+e", this->uene, this->vary, this->alph); // 14
        fprintf(fp, " %+e %+e %+e", this->dens, this->vsnd, this->pres); // 17
        fprintf(fp, " %+e %+e %+e", this->divv, this->rotv, this->bswt); // 20
        fprintf(fp, " %+e %+e %+e", this->grdh, this->ksr,  this->pot); // 23
        fprintf(fp, " %6d %6d", this->np, this->np2); // 25
        fprintf(fp, " %+e %+e %+e", this->vsmx, this->udot, this->ydot);
        fprintf(fp, "\n");
    }

    void referEquationOfState() {
        this->pres  = (RP::AdiabaticIndex - 1.) * this->dens * this->uene;
        this->presk = pow(this->pres, RP::PowerForWeight);
        this->vsnd  = sqrt(RP::AdiabaticIndex * this->pres / this->dens);
        this->dpdr  = (RP::AdiabaticIndex - 1.) * this->uene;
        this->dpdu  = (RP::AdiabaticIndex - 1.) * this->dens;
        this->gamm  = RP::AdiabaticIndex;
    }

    PS::F64 calcTimeStep() {
        PS::F64 tceff   = RP::CoefficientOfTimestep;
        PS::F64 dth = tceff * this->ksr / (this->vsmx * SK::ksrh);
        PS::F64 dtu = tceff * fabs(this->uene / this->udot);
        PS::F64 dty = tceff * fabs(this->vary / this->ydot);
        return std::min(std::min(dth, dtu), dty);
    }

    PS::F64 calcEnergy() {
        return this->mass * (0.5 * this->vel * this->vel + this->uene + 0.5 * this->pot);
    }

    PS::F64vec calcMomentum() {
        return this->mass * this->vel;
    }

    void predict(PS::F64 dt) {
        PS::F64 pdot = - this->dpdr * this->dens * this->divv + this->dpdu * this->udot;
        PS::F64 hdot =   this->ksr  * this->divv / RP::NumberOfDimension;
        this->pos    = this->pos   +       this->vel   * dt  + 0.5 * this->acc * dt * dt;
        this->vel2   = this->vel   + 0.5 * this->acc   * dt;
        this->vel    = this->vel   +       this->acc   * dt;
        this->uene2  = this->uene  + 0.5 * this->udot  * dt;
        this->uene   = this->uene  +       this->udot  * dt;
        //this->vary2  = this->vary  + 0.5 * this->ydot  * dt;
        this->vary   = this->vary  +       this->ydot  * dt;
        this->alph2  = this->alph  + 0.5 * this->adot  * dt;
        this->alph   = this->alph  +       this->adot  * dt;
        this->pres   = this->pres  +             pdot  * dt;
        this->ksr    = this->ksr   +             hdot  * dt;
        this->presk  = pow(this->pres, RP::PowerForWeight);
    }

    void correct(PS::F64 dt) {
        this->vel   = this->vel2   + 0.5 * this->acc   * dt;
        this->uene  = this->uene2  + 0.5 * this->udot  * dt;
        //this->vary  = this->vary2  + 0.5 * this->ydot  * dt;
        this->alph  = this->alph2  + 0.5 * this->adot  * dt;
    }
};

namespace RunParameter {
    PS::S64 readAscii(const char * const filename) {
        FILE *fp = fopen(filename, "r");
        fscanf(fp, "%lf%lf%lf", &Time, &TimeEnd, &TimestepAscii);
        fscanf(fp, "%lld%lld", &NumberOfDimension, &KernelType);
        fscanf(fp, "%lf%lf%lf", &ComputationalBox.low_[0], &ComputationalBox.low_[1], &ComputationalBox.low_[2]);
        fscanf(fp, "%lf%lf%lf", &ComputationalBox.high_[0], &ComputationalBox.high_[1], &ComputationalBox.high_[2]);
        fscanf(fp, "%lf%lf%lf%lf", &AdiabaticIndex, &AlphaMaximum, &AlphaMinimum, &CoefficientOfTimestep);
        fscanf(fp, "%lld", &NumberOfParticle);
        fclose(fp);
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
          class Tvolume,
          class Tauxiliary,
          class Thydro>
void startSimulation(char **argv,
                     Tdinfo & dinfo,
                     Tsph & sph,
                     Tvolume & volume,
                     Tauxiliary & auxiliary,
                     Thydro & hydro) {
    RP::readAscii(argv[2]);
    RP::KernelSupportRadiusMaximum = std::numeric_limits<PS::F64>::max();
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

    PS::TreeForForceShort<Density, DensityEPI, DensityEPJ>::Gather density;
    density.initialize(0);
    calcDensityKernel(dinfo, sph, density);
    referEquationOfState(sph);
    calcVariableY(sph);
    calcSPHKernel(dinfo, sph, volume, auxiliary, hydro);
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
          class Tvolume,
          class Tauxiliary,
          class Thydro>
void loopSimulation(Tdinfo & dinfo,
                    Tsph & sph,
                    Tvolume & volume,
                    Tauxiliary & auxiliary,
                    Thydro & hydro) {

    WT::clear();
    while(RP::Time < RP::TimeEnd) {
        RP::NumberOfStep++;
        /*
        if(RP::Time > 0.) {
            sph.writeParticleAscii("snap/hoge.dat");
            PS::Finalize();
            exit(0);
        }
        */
        RP::Timestep = calcTimeStep(sph);
        WT::reduceInterProcess();
        outputData(sph);
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
        calcSPHKernel(dinfo, sph, volume, auxiliary, hydro);
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
