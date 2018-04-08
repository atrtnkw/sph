#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>
#include <mpi.h>

//#define DEBUG

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};

#include "particle_simulator.hpp"
#include "hdr_time.hpp"
#include "hdr_run.hpp"
#ifdef USE_INTRINSICS
#include "vector_x86.hpp"
#endif
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_hgas.hpp"

class SPH3D : public HelmholtzGas {
public:
    SPH3D() {
        this->id    = 0;
        this->istar = 0;
        this->mass  = 0.;
        this->pos   = 0.;
        this->vel   = 0.;
        this->uene  = 0.;
        this->alph  = 0.;
        this->alphu = 0.;
        this->ksr   = 0.;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            this->cmps[k] = 0.;
        }        
    }

    void readAscii(FILE * fp) {
        fscanf(fp, "%lld%lld%lf", &this->id, &this->istar, &this->mass);      //  3
        fscanf(fp, "%lf%lf%lf", &this->pos[0], &this->pos[1], &this->pos[2]); //  6
        fscanf(fp, "%lf%lf%lf", &this->vel[0], &this->vel[1], &this->vel[2]); //  9
        fscanf(fp, "%lf%lf%lf", &this->acc[0], &this->acc[1], &this->acc[2]); // 12
        fscanf(fp, "%lf%lf%lf", &this->uene, &this->alph, &this->alphu);      // 15
        fscanf(fp, "%lf%lf%6d", &this->dens, &this->ksr,  &this->np);         // 18
        fscanf(fp, "%lf%lf%lf", &this->vsnd, &this->pres, &this->temp);       // 21
        fscanf(fp, "%lf%lf%lf", &this->divv, &this->rotv, &this->bswt);       // 24
        fscanf(fp, "%lf%lf%lf", &this->pot,  &this->abar, &this->zbar);       // 27
        fscanf(fp, "%lf",       &this->enuc);                                 // 28
        fscanf(fp, "%lf%lf%lf", &this->vsmx, &this->udot, &this->dnuc);       // 31
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) { // 32 -- 44
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        fscanf(fp, "%lf", &this->pot3);
        fscanf(fp, "%lf%lf%lf", &this->tempmax[0], &this->tempmax[1], &this->tempmax[2]);
        fscanf(fp, "%lf", &this->entr);
    }

    void writeAscii(FILE * fp = stdout) const {
        fprintf(fp, "%8d %2d %+e",  this->id, this->istar, this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->uene, this->alph, this->alphu);
        fprintf(fp, " %+e", this->ksr);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
#ifdef DEBUG
        fprintf(fp, " %+e %+e %+e", this->dens, this->temp);
#endif
        fprintf(fp, "\n");
    }

    void setNonZeroToNucleus() {
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            this->cmps[k] = ((this->cmps[k] != 0.) ? this->cmps[k] : 1e-10);
        }
    }

};

namespace HeliumDetonation {
    PS::F64    tmax  = 1.0e9;
    PS::F64    tmin  = 1.0e7;
    PS::F64vec nvec  = 0.;
    PS::F64    width = 0.;
    bool       first = true;

    PS::F64 getTemperatureLinear(SPH3D & sph) {
        PS::F64 temp = 0.;
        
        if(first) {
            PS::F64 linv = 1. / sqrt(nvec * nvec);
            nvec = nvec * linv;
            first = false;
        }

        PS::F64 dist = nvec * sph.pos;
        if(dist > width) {
            temp = tmax;
        } else if (dist > 0.) {
            temp = tmin + (tmax - tmin) * dist / width;
        } else {
            temp = tmin;
        }

        return temp;
    }
}

int main(int argc, char ** argv)
{
    PS::Initialize(argc, argv);
    
    init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);

    PS::ParticleSystem<SPH3D> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);

    char itype[1024];
    char otype[1024];
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s%s", itype, otype);
    {
        using namespace HeliumDetonation;
        fscanf(fp, "%lf%lf%lf%lf", &nvec[0], &nvec[1], &nvec[2], &width);
    }
    if(PS::Comm::getRank() == 0) {
        printf("itype: %s\n", itype);
        printf("otype: %s\n", otype);
    }
    fclose(fp);

#ifdef DEBUG
    char dzfile[1024];
    sprintf(dzfile, "%s_p%06d_i%06d.dat.dz", otype, PS::Comm::getNumberOfProc(),
            PS::Comm::getRank());
    FILE * zfp = fopen(dzfile, "w");
#endif

    sph.readParticleAscii(itype, "%s_p%06d_i%06d.dat");
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].temp = HeliumDetonation::getTemperatureLinear(sph[i]);
        sph[i].setNonZeroToNucleus();
        flash_helmholtz_e_(&sph[i].dens, &sph[i].temp, sph[i].cmps.getPointer(), &sph[i].uene);

#ifdef DEBUG
        if(-1e4 < sph[i].pos[2] && sph[i].pos[2] < 1e4) {
            sph[i].writeAscii(zfp);
        }
#endif
    }
#ifdef DEBUG
    fclose(zfp);
#endif
    sph.writeParticleAscii(otype, "%s_p%06d_i%06d.data");


    PS::Finalize();

    return 0;
}
