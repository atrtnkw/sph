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
    const  PS::S64    NumberOfMesh = 256;
    static PS::F64    dens_l[NumberOfMesh], dens_g[NumberOfMesh];
    static PS::F64    temp_l[NumberOfMesh], temp_g[NumberOfMesh];
    static PS::F64    vsnd_l[NumberOfMesh], vsnd_g[NumberOfMesh];
    static PS::F64    vphi_l[NumberOfMesh], vphi_g[NumberOfMesh];
    //static PS::F64    cmps_l[NumberOfMesh][NR::NumberOfNucleon], cmps_g[NumberOfMesh][NR::NumberOfNucleon];
    static PS::F64    rphi[NumberOfMesh];
    static PS::F64vec mesh[NumberOfMesh];

    PS::S64 TargetId;
    PS::F64 TargetRadius;
    PS::F64 Phi0;
    PS::F64 Phi1;

    PS::F64 CoreRadius;
    PS::F64 ShellThickness;
    PS::F64 DeltaTheta;
    PS::F64 MinimumPhi;
    PS::F64 MaximumPhi;
    PS::F64 CarbonFraction;

    void putMeshPoint() {
        PS::F64 dPhi = (Phi1 - Phi0) / (PS::F64)NumberOfMesh;
        for(PS::S64 i = 0; i < NumberOfMesh; i++) {
            PS::F64 iPhi   = Phi0 + dPhi * (i + 0.5);
            rphi[i]    = TargetRadius * iPhi;
            mesh[i][0] = TargetRadius * cos(iPhi);
            mesh[i][1] = TargetRadius * sin(iPhi);
            mesh[i][2] = 0.;
        }
    }

    void readRunParameter(FILE * fp) {
        fscanf(fp, "%lld", &TargetId);
        fscanf(fp, "%lf", &TargetRadius);
        fscanf(fp, "%lf%lf", &Phi0, &Phi1);
        Phi0 *= M_PI / 180.;
        Phi1 *= M_PI / 180.;
        putMeshPoint();
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
        fscanf(fp, "%lf%lf%6d", &this->dens, &this->ksr,  &this->np);         // 18
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
void obtainDensityCenter(Tsph & sph,
                         PS::F64vec & pcenter,
                         PS::F64vec & vcenter) {
    using namespace TargetWhiteDwarf;

    PS::S64 NotTargetId = (TargetId + 1) % 2;

    PS::S64    np   = 0;
    PS::F64    dloc = 0.;
    PS::F64vec ploc(0.);
    PS::F64vec vloc(0.);
    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        if(sph[ip].istar == NotTargetId) { // Temporary
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
void extractHeDetonation(char * ofile,
                         Tsph & sph) {
    using namespace TargetWhiteDwarf;

    for(PS::S64 i = 0; i < NumberOfMesh; i++) {
        dens_l[i] = 0.;
        temp_l[i] = 0.;
        vsnd_l[i] = 0.;
        vphi_l[i] = 0.;
    }

    PS::F64vec pcenter(0.);
    PS::F64vec vcenter(0.);
    obtainDensityCenter(sph, pcenter, vcenter);

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

            PS::F64 ivphi = (sph[ip].pos[0] * sph[ip].vel[1] - sph[ip].pos[1] * sph[ip].vel[0])
                / sqrt(sph[ip].pos[0] * sph[ip].pos[0] + sph[ip].pos[1] * sph[ip].pos[1]);
            vphi_l[jp] += ksph * ivphi;
        }

    } 

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(dens_l, dens_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(temp_l, temp_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(vsnd_l, vsnd_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(vphi_l, vphi_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        for(PS::S64 i = 0; i < NumberOfMesh; i++) {
            PS::F64 xlen = rphi[i] - TargetRadius * Phi0;
            fprintf(fp, " %+e %+e %+e", xlen,      mesh[i][0], mesh[i][1]);
            fprintf(fp, " %+e %+e %+e", dens_g[i], temp_g[i],  vsnd_g[i]);
            fprintf(fp, " %+e"        , vphi_g[i]);
            fprintf(fp, "\n");
        }
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
    PS::S64 ibgn, iend;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld%lld", &ibgn, &iend);
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
        extractHeDetonation(ofile, sph);

    }

    PS::Finalize();

    return 0;
}
