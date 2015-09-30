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
        fprintf(fp, "%6d %+e", this->id, this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->acc[0], this->acc[1], this->acc[2]);
        fprintf(fp, " %+e %+e %+e", this->uene, this->alph, this->ksr);
        fprintf(fp, " %+e %+e %+e", this->dens, this->vsnd, this->pres);
        fprintf(fp, " %+e %+e %+e", this->divv, this->rotv, this->bswt);
        fprintf(fp, " %+e %6d %+e", this->grdh, this->np, this->pot);
        fprintf(fp, " %+e", this->eta);
        ///////////////////
        fprintf(fp, " %+e %+e %+e", this->ctau.xx, this->ctau.yy, this->ctau.zz);
        fprintf(fp, " %+e %+e %+e", this->ctau.xy, this->ctau.xz, this->ctau.yz);
        ///////////////////
        fprintf(fp, "\n");
    }

    void referEquationOfState() {
        this->pres = (RP::AdiabaticIndex - 1.)  * this->dens * this->uene;
        this->vsnd = sqrt(RP::AdiabaticIndex * this->pres / this->dens);
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

template <class Tdinfo,
          class Tsph,
          class Tvolume,
          class Tauxiliary>
void startSimulation(char **argv,
                     Tdinfo & dinfo,
                     Tsph & sph,
                     Tvolume & volume,
                     Tauxiliary & auxiliary) {
    RP::readAscii(argv[2]);
    RP::KernelSupportRadiusMaximum = std::numeric_limits<PS::F64>::max();
    sph.readParticleAscii(argv[3]);    
    if(RP::ComputationalBox.low_[0] != RP::ComputationalBox.high_[0]) {
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
        dinfo.setPosRootDomain(RP::ComputationalBox.low_, RP::ComputationalBox.high_);
    }
    ND::setDimension(RP::NumberOfDimension);
    SK::setKernel(RP::KernelType, RP::NumberOfDimension);
    WT::start();
    sph.adjustPositionIntoRootDomain(dinfo);
    WT::accumulateOthers();
    WT::start();
    PS::MT::init_genrand(0);
    dinfo.decomposeDomainAll(sph);
    WT::accumulateDecomposeDomain();
    WT::start();
    sph.exchangeParticle(dinfo);
    WT::accumulateExchangeParticle();

    PS::TreeForForceShort<Density, DensityEPI, DensityEPJ>::Gather density;
    density.initialize(0);
    calcDensityKernel(dinfo, sph, density);
    referEquationOfState(sph);
    calcSPHKernel(dinfo, sph, volume, auxiliary);
}

template <class Tdinfo,
          class Tsph>
void restartSimulation(char **argv,
                       Tdinfo & dinfo,
                       Tsph & sph) {
    return;
}

typedef IdealGas GeneralSPH;
