#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>
#include <mpi.h>

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
#include "hdr_bhns.hpp"

//MPI_Datatype NucleonMPI;

static PS::S64 FlagOfBoundOrNot = -2;
// -1: All particles
//  0: Bound particles
//  1: Unbound particles
static PS::S64 FlagOfParticle   = -2;
// -1: All particles
//  0: Primary particle
//  1: Secondary particle

class SPHAnalysis : public HelmholtzGas {
public:
    SPHAnalysis() {
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
        fscanf(fp, "%lf%lf%lld", &this->dens, &this->ksr,  &this->np);        // 18
        fscanf(fp, "%lf%lf%lf", &this->vsnd, &this->pres, &this->temp);       // 21
        fscanf(fp, "%lf%lf%lf", &this->divv, &this->rotv, &this->bswt);       // 24
        fscanf(fp, "%lf%lf%lf", &this->pot,  &this->abar, &this->zbar);       // 27
        fscanf(fp, "%lf",       &this->enuc);                                 // 28
        fscanf(fp, "%lf%lf%lf", &this->vsmx, &this->udot, &this->dnuc);       // 31
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) {       // 32 -- 44
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        fscanf(fp, "%lf", &this->pot3);
        fscanf(fp, "%lf%lf%lf", &this->tempmax[0], &this->tempmax[1], &this->tempmax[2]);
        fscanf(fp, "%lf", &this->entr);
    }

};

template <class Tsph>
void obtainDensityCenter(Tsph & sph,
                         PS::F64vec & pcenter,
                         PS::F64vec & vcenter) {
    PS::S64    np   = 0;
    PS::F64    dloc = 0.;
    PS::F64vec ploc(0.);
    PS::F64vec vloc(0.);
    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        if(sph[ip].istar == 0) {
            continue;
        }
        np++;
        dloc += sph[ip].dens;
        ploc += sph[ip].dens * sph[ip].pos;
        vloc += sph[ip].dens * sph[ip].vel;
    }

    PS::F64    dglb;
    PS::F64vec pglb;
    PS::F64vec vglb;
    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(&dloc, &dglb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&ploc, &pglb, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&vloc, &vglb, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    pcenter = pglb / dglb;
    vcenter = vglb / dglb;
}

template <class Tsph>
void reduceElement(Tsph & sph,
                   PS::F64 * cmpsglb) {

    PS::F64vec pcenter(0.);
    PS::F64vec vcenter(0.);
    obtainDensityCenter(sph, pcenter, vcenter);

    PS::F64 cmpsloc[NR::NumberOfNucleon];
    for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
        cmpsloc[k] = 0.;
        cmpsglb[k] = 0.;
    }
    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        PS::F64vec pptcl = sph[ip].pos - pcenter;
        PS::F64vec vptcl = sph[ip].vel - vcenter;
        PS::F64    eptcl = 0.5 * (vptcl * vptcl) + sph[ip].pot;

        if((FlagOfBoundOrNot == 0 && eptcl > 0.)
           || (FlagOfBoundOrNot == 1 && eptcl <= 0.)) {
            continue;
        }

        if((FlagOfParticle != -1) && (FlagOfParticle != sph[ip].istar)) {
            continue;
        }

        for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
            cmpsloc[k] += sph[ip].cmps[k] * sph[ip].mass;
        }
    }

    PS::S64 ret = MPI_Allreduce(cmpsloc, cmpsglb, NR::NumberOfNucleon,
                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    assert(ret == MPI_SUCCESS);

    for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
        cmpsglb[k] = cmpsglb[k] / CodeUnit::SolarMass;
    }

    return;
}

int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);

    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);

    char idir[1024], odir[1024];
    PS::S64 ibgn, iend;
    PS::F64 tsnap;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld%lld%lf", &ibgn, &iend, &tsnap);
    fscanf(fp, "%lld", &FlagOfBoundOrNot);
    fscanf(fp, "%lld", &FlagOfParticle);
    fclose(fp);

    assert(-1 <= FlagOfBoundOrNot && FlagOfBoundOrNot <= 1);
    assert(-1 <= FlagOfParticle && FlagOfParticle <= 1);

    char ofile[1024];
    sprintf(ofile, "%s/element.dat", odir);
    FILE * fpout = fopen(ofile, "w");

    for(PS::S64 itime = ibgn; itime <= iend; itime++) {        
        char tfile[1024];
        FILE *fp = NULL;
        PS::S64 tdir = 0;
        for(PS::S64 iidir = 0; iidir < 100; iidir++) {
            sprintf(tfile, "%s/t%02d/sph_t%04d_p%06d_i%06d.dat", idir, iidir, itime,
                    PS::Comm::getNumberOfProc(), 0);
            fp = fopen(tfile, "r");
            if(fp != NULL) {
                tdir = iidir;
                break;
            }
        }
        if(fp == NULL) {
            if(PS::Comm::getRank() == 0) {
                fprintf(stderr, "Not found %s\n", tfile);
            }
            continue;
        }
        fclose(fp);

        char sfile[1024];
        sprintf(sfile, "%s/t%02d/sph_t%04d", idir, tdir, itime);
        sph.readParticleAscii(sfile, "%s_p%06d_i%06d.dat");

        PS::F64 cmps[NR::NumberOfNucleon];
        reduceElement(sph, cmps);
        if(PS::Comm::getRank() == 0) {
            fprintf(fpout, "%+e", itime * tsnap);
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                fprintf(fpout, " %+e", cmps[k]);
            }
            fprintf(fpout, "\n");
            fflush(fpout);
        }

    }

    fclose(fpout);

    PS::Finalize();

    return 0;
}
