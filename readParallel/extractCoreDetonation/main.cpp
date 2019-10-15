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

namespace TargetWhiteDwarf{
    const  PS::S64    NumberOfMesh = 1024;
    static PS::F64    dens_l[NumberOfMesh], dens_g[NumberOfMesh];
    static PS::F64    temp_l[NumberOfMesh], temp_g[NumberOfMesh];
    static PS::F64    vsnd_l[NumberOfMesh], vsnd_g[NumberOfMesh];
    static PS::F64    vrad_l[NumberOfMesh], vrad_g[NumberOfMesh];
    static PS::F64    cmps_l[NumberOfMesh][NR::NumberOfNucleon];
    static PS::F64    cmps_g[NumberOfMesh][NR::NumberOfNucleon];
    static PS::F64vec mesh[NumberOfMesh];

    PS::S64 TargetId;
    PS::F64 Radius;
    PS::F64 Phi;
    PS::F64 M0, M1;
    PS::F64 BinaryRadius;
    PS::F64 BinaryVelocity;
    PS::F64 BinaryOmega;

    void putMeshPoint() {
#if 1
        PS::F64 dRadius = (2. * Radius) / NumberOfMesh;
        for(PS::S64 i = 0; i < NumberOfMesh; i++) {
            mesh[i][0] = cos(Phi) * ((0.5 + (PS::F64)i) * dRadius - Radius);
            mesh[i][1] = sin(Phi) * ((0.5 + (PS::F64)i) * dRadius - Radius);
            mesh[i][2] = 0.;
        }
#else
        PS::F64 dRadius = Radius / NumberOfMesh;
        for(PS::S64 i = 0; i < NumberOfMesh; i++) {
            mesh[i][0] = cos(Phi) * ((0.5 + (PS::F64)i) * dRadius);
            mesh[i][1] = sin(Phi) * ((0.5 + (PS::F64)i) * dRadius);
            mesh[i][2] = 0.;
        }
#endif
    }

    void readRunParameter(FILE * fp) {
        fscanf(fp, "%lld", &TargetId);
        fscanf(fp, "%lf%lf", &Radius, &Phi);
        fscanf(fp, "%lf%lf", &M0, &M1);
        fscanf(fp, "%lf", &BinaryRadius);
        fscanf(fp, "%lf", &BinaryVelocity);
        Phi        *= M_PI / 180.;
        BinaryOmega = BinaryVelocity / BinaryRadius;
        putMeshPoint();
    }

    void obtainCoordinateCenter(PS::F64 time,
                                PS::F64vec & pcenter,
                                PS::F64vec & vcenter) {
        if(TargetId == 0) {
            PS::F64 ratio = M1 / (M0 + M1);
            pcenter[0] = + ratio * BinaryRadius * cos(BinaryOmega * time);
            pcenter[1] = + ratio * BinaryRadius * sin(BinaryOmega * time);
            pcenter[2] = 0.;
            vcenter[0] = - ratio * BinaryRadius * sin(BinaryOmega * time) * BinaryOmega;
            vcenter[1] = + ratio * BinaryRadius * cos(BinaryOmega * time) * BinaryOmega;
            vcenter[2] = 0.;
        } else {
            PS::F64 ratio = M0 / (M0 + M1);
            pcenter[0] = - ratio * BinaryRadius * cos(BinaryOmega * time);
            pcenter[1] = - ratio * BinaryRadius * sin(BinaryOmega * time);
            pcenter[2] = 0.;
            vcenter[0] = + ratio * BinaryRadius * sin(BinaryOmega * time) * BinaryOmega;
            vcenter[1] = - ratio * BinaryRadius * cos(BinaryOmega * time) * BinaryOmega;
            vcenter[2] = 0.;
        }
    }

}

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
//        fscanf(fp, "%lf%lf%6d", &this->dens, &this->ksr,  &this->np);         // 18
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

    void writeAscii(FILE * fp) {
        fprintf(fp, "%6d %2d %+e",  this->id,     this->istar,  this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->dens,   this->temp,   this->vsnd);
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) {
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        fprintf(fp, "\n");
    }

};

template <class Tsph>
void extractHeDetonation(char * ofile,
                         PS::F64 time,
                         Tsph & sph) {
    using namespace TargetWhiteDwarf;

    for(PS::S64 i = 0; i < NumberOfMesh; i++) {
        dens_l[i] = 0.;
        temp_l[i] = 0.;
        vsnd_l[i] = 0.;
        vrad_l[i] = 0.;
        for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
            cmps_l[i][k] = 0.;
        }
    }

    PS::F64vec pcenter(0.);
    PS::F64vec vcenter(0.);
    obtainCoordinateCenter(time, pcenter, vcenter);

    if(PS::Comm::getRank() == 0) {
        printf("t: %+e\n", time);
        printf("x: %+e %+e %+e\n", pcenter[0], pcenter[1], pcenter[2]);
        printf("v: %+e %+e %+e\n", vcenter[0], vcenter[1], vcenter[2]);
    }

    PS::F64vec raxis(cos(TargetWhiteDwarf::Phi), sin(TargetWhiteDwarf::Phi), 0.);
    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        sph[ip].pos -= pcenter;
        sph[ip].vel -= vcenter;

        for(PS::S64 jp = 0; jp < NumberOfMesh; jp++) {
            PS::F64vec dr    = mesh[jp] - sph[ip].pos;
            PS::F64    dr2   = dr * dr;
            PS::F64    dr1   = sqrt(dr2);
            PS::F64    hinv1 = 1. / sph[ip].ksr;
            PS::F64    hinv3 = ND::calcVolumeInverse(hinv1);
            PS::F64    dq1   = dr1 * hinv1;
            PS::F64    kw0   = hinv3 * SK::kernel0th(dq1);
            dens_l[jp] += sph[ip].mass * kw0;

            PS::F64 dinv = 1. / sph[ip].dens;
            PS::F64 ksph = sph[ip].mass * dinv * kw0;
            temp_l[jp] += ksph * sph[ip].temp;
            vsnd_l[jp] += ksph * sph[ip].vsnd;
            vrad_l[jp] += ksph * (raxis * sph[ip].vel);
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                cmps_l[jp][k] += ksph * sph[ip].cmps[k];
            }
        }

    } 

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(dens_l, dens_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(temp_l, temp_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(vsnd_l, vsnd_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(vrad_l, vrad_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(cmps_l, cmps_g, NumberOfMesh*NR::NumberOfNucleon, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
#if 1
        PS::F64 dRadius = (2. * Radius) / NumberOfMesh;
        for(PS::S64 i = 0; i < NumberOfMesh; i++) {
            PS::F64 xlen = (0.5 + (PS::F64)i) * dRadius - Radius;
            PS::F64 cmps = 0.;
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                cmps += cmps_g[i][k];
            }
            fprintf(fp, " %+e %+e %+e", xlen,      mesh[i][0], mesh[i][1]);
            fprintf(fp, " %+e %+e %+e", dens_g[i], temp_g[i],  vsnd_g[i]);
            fprintf(fp, " %+e %+e %+e"    , vrad_g[i], cmps_g[i][0]/cmps, cmps_g[i][1]/cmps);
            fprintf(fp, "\n");
        }
#else
        PS::F64 dRadius = Radius / NumberOfMesh;
        for(PS::S64 i = 0; i < NumberOfMesh; i++) {
            PS::F64 xlen = (0.5 + (PS::F64)i) * dRadius;
            fprintf(fp, " %+e %+e %+e", xlen,      mesh[i][0], mesh[i][1]);
            fprintf(fp, " %+e %+e %+e", dens_g[i], temp_g[i],  vsnd_g[i]);
            fprintf(fp, " %+e"        , vrad_g[i]);
            fprintf(fp, "\n");
        }
#endif
        fclose(fp);
    }

}

int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);

    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);

    char idir[1024], odir[1024];
    PS::S64 ibgn, iend, nsnap;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld%lld%lld", &ibgn, &iend, &nsnap);
    TargetWhiteDwarf::readRunParameter(fp);
    fclose(fp);

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

        char ofile[1024];
        sprintf(ofile, "%s/mesh_t%04d.dat", odir, itime);
        extractHeDetonation(ofile, itime/(PS::F64)nsnap, sph);

    }

    PS::Finalize();

    return 0;
}
