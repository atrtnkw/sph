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

namespace AssignToMesh {
    const PS::S64 NumberOfMesh    = 800;
    const PS::F64 HalfLengthOfBox = 5e8;
    const PS::F64 WidthOfMesh     = (2. * HalfLengthOfBox) / ((PS::F64)NumberOfMesh);
    const PS::F64 WidthOfMeshInv  = 1. / WidthOfMesh;
    static PS::F64 dens_l[NumberOfMesh];
    static PS::F64 dens_g[NumberOfMesh];
    static PS::F64 velz_l[NumberOfMesh];
    static PS::F64 velz_g[NumberOfMesh];

    PS::S64 getIndex(PS::F64 x) {
        PS::F64 fx = (x + HalfLengthOfBox) * WidthOfMeshInv - 0.5;
        PS::S64 ix = (PS::S64)(ceil(fx)); // round-up
        ix = (ix >= 0)           ? ix : 0;
        ix = (ix < NumberOfMesh) ? ix : (NumberOfMesh - 1);
        return ix;
    }

    PS::F64 getPosition(PS::S64 index) {
        return (- HalfLengthOfBox + WidthOfMesh * ((PS::F64)index + 0.5));
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
void assignToMesh(char * ofile,
                  Tsph & sph) {
    using namespace AssignToMesh;

    for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
        dens_l[iz] = 0.;
        velz_l[iz] = 0.;
    }

    PS::F64 posx = +5.854266e+08;
    PS::F64 posy = -1.250326e+09;
//    PS::F64 posx = +1.594931e+08;
//    PS::F64 posy = -1.375774e+09;
    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
            PS::F64vec pos(posx, posy, getPosition(iz));
            PS::F64vec dr = pos - sph[ip].pos;
            PS::F64 dr2   = dr * dr;
            PS::F64 dr1   = sqrt(dr2);
            PS::F64 hinv1 = 1. / sph[ip].ksr;
            PS::F64 hinv3 = ND::calcVolumeInverse(hinv1);
            PS::F64 dq1   = dr1 * hinv1;
            PS::F64 kw0   = hinv3 * SK::kernel0th(dq1);
            dens_l[iz] += sph[ip].mass * kw0;
            PS::F64 dinv = 1. / sph[ip].dens;
            PS::F64 ksph = sph[ip].mass * dinv * kw0;
            velz_l[iz] += ksph * sph[ip].vel[2];
        }
    }

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(dens_l, dens_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(velz_l, velz_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
            PS::F64 mz   = getPosition(iz);
            PS::F64 dens = (dens_g[iz] > 1e-6) ? dens_g[iz] : 1e-6;
            fprintf(fp, "%+e %+e %+e\n", mz, dens, velz_g[iz]);
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

    char idir[1024], odir[1024], otype[1024];
    PS::S64 ibgn, iend;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%s", otype);
    fscanf(fp, "%lld%lld", &ibgn, &iend);
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
        sprintf(ofile, "%s/%s.dat", odir, otype);
        assignToMesh(ofile, sph);

    }

    PS::Finalize();

    return 0;
}
