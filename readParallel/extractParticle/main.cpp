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

MPI_Datatype NucleonMPI;

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

int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);

    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);

    char idir[1024], odir[1024];
    PS::S64 ibgn, iend;
    PS::S64 itime;
    PS::S64 id1st, idint;
    PS::S64 mode;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld", &itime);
    fscanf(fp, "%lld%lld", &id1st, &idint);
    fscanf(fp, "%lld", &mode);
    fclose(fp);

    char tfile[1024];
    fp = NULL;
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
        PS::Finalize();
        exit(0);
    }
    fclose(fp);
    
    char sfile[1024];
    sprintf(sfile, "%s/t%02d/sph_t%04d", idir, tdir, itime);
    sph.readParticleAscii(sfile, "%s_p%06d_i%06d.dat");
    
    char ofile[1024];
    sprintf(ofile, "%s/ext_t%04d_p%06d_i%06d.dat", odir, itime, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
    fp = fopen(ofile, "w");
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        if((sph[i].id-id1st)%idint != 0) {
            continue;
        }
        if(mode == 0) {
            fprintf(fp, "%10d %2d %+e", sph[i].id, sph[i].istar, sph[i].mass);
            fprintf(fp, " %+e %+e %+e", sph[i].pos[0], sph[i].pos[1], sph[i].pos[2]);
            fprintf(fp, " %+e %+e %+e", sph[i].vel[0], sph[i].vel[1], sph[i].vel[2]);
            fprintf(fp, " %+e %+e %+e", sph[i].dens, sph[i].ksr, sph[i].pot);
            fprintf(fp, "\n");
        } else if(mode == 1) {
            if(sph[i].istar == 0 && sph[i].cmps[0] > 0.1) {
                fprintf(fp, "%10d\n", sph[i].id);
            }
        } else {
            if(sph[i].istar == 0) {
                sph[i].writeAscii(fp);
            }
        }
    }
    fclose(fp);

    PS::Finalize();

    return 0;
}
