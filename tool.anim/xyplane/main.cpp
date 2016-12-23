#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>
#include <mpi.h>
//#include <Eigen/Eigenvalues>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};

#include "particle_simulator.hpp"
#include "hdr_time.hpp"
#include "hdr_run.hpp"
#include "vector_x86.hpp"
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_hgas.hpp"
#include "hdr_bhns.hpp"

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

    /*
    void write(FILE * fp = stdout) {
        fprintf(fp, "%8d %2d %+e",  this->id, this->istar, this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->uene, this->alph, this->alphu);
        fprintf(fp, " %+e", this->ksr);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        fprintf(fp, "\n");
    }
    */

    inline PS::F64mat calcMomentOfInertia(PS::F64vec xc) {
        PS::F64mat qq = 0.;
        PS::F64vec dx = this->pos - xc;
        PS::F64 r2 = dx * dx;
        qq.xx = this->mass * (r2 - dx[0] * dx[0]);
        qq.xy = this->mass * (   - dx[0] * dx[1]);
        qq.xz = this->mass * (   - dx[0] * dx[2]);
        qq.yy = this->mass * (r2 - dx[1] * dx[1]);
        qq.yz = this->mass * (   - dx[1] * dx[2]);
        qq.zz = this->mass * (r2 - dx[2] * dx[2]);
        return qq;
    }
};


class BlackHoleNeutronStarAnalysis : public BlackHoleNeutronStar {
public:
    BlackHoleNeutronStarAnalysis() {
        this->id    = 0;
        this->istar = 0;
        this->mass  = 0.;
        this->pos   = 0.;
        this->vel   = 0.;
        this->acc   = 0.;
        this->eps   = 0.;
        this->pot   = 0.;
    }

    void readAscii(FILE * fp) {
        fscanf(fp, "%lld%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->istar, &this->mass,
               &this->pos[0], &this->pos[1], &this->pos[2],
               &this->vel[0], &this->vel[1], &this->vel[2],
               &this->acc[0], &this->acc[1], &this->acc[2],
               &this->eps, &this->pot);
    }

    inline PS::F64mat calcMomentOfInertia(PS::F64vec xc) {
        PS::F64mat qq = 0.;
        PS::F64vec dx = this->pos - xc;
        PS::F64 r2 = dx * dx;
        qq.xx = this->mass * (r2 - dx[0] * dx[0]);
        qq.xy = this->mass * (   - dx[0] * dx[1]);
        qq.xz = this->mass * (   - dx[0] * dx[2]);
        qq.yy = this->mass * (r2 - dx[1] * dx[1]);
        qq.yz = this->mass * (   - dx[1] * dx[2]);
        qq.zz = this->mass * (r2 - dx[2] * dx[2]);
        return qq;
    }
};

template <class Tsph,
          class Tbhns>
void projectOnOrbitalPlane(char * ofile,
                           Tsph  & sph,
                           Tbhns & bhns,
                           PS::F64vec xmin,
                           PS::F64    wdth,
                           PS::S64    nmax) {

    PS::F64 dx    = wdth / (PS::F64)nmax;
    PS::F64 dxinv = 1. / dx;

    PS::S64 ptcl[nmax][nmax];
    PS::F64 dens[nmax][nmax];
    PS::F64 temp[nmax][nmax];
    PS::F64 entr[nmax][nmax];

    for(PS::S64 i = 0; i < nmax; i++) {
        for(PS::S64 j = 0; j < nmax; j++) {
            ptcl[i][j] = 0;
            dens[i][j] = 0.;
            temp[i][j] = 0.;
            entr[i][j] = 0.;
        }
    }

    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
#if 0
        if(fabs(sph[i].pos[2]) > fabs(sph[i].ksr)) {
            continue;
        }
        PS::F64vec px = sph[i].pos - bhns[0].pos;
#else
        PS::F64vec px = sph[i].pos - bhns[0].pos;
        if(fabs(px[2]) > fabs(sph[i].ksr)) {
            continue;
        }
#endif
        PS::S64 nxi = (PS::S64)((px[0] - xmin[0]) * dxinv);
        if(nxi < 0 || nmax <= nxi) {
            continue;
        }
        PS::S64 nyi = (PS::S64)((px[1] - xmin[1]) * dxinv);
        if(nyi < 0 || nmax <= nyi) {
            continue;
        }
        ptcl[nxi][nyi] += 1;
        dens[nxi][nyi] += sph[i].dens;
        temp[nxi][nyi] += sph[i].temp;
        entr[nxi][nyi] += sph[i].entr;
    }

    for(PS::S64 i = 0; i < nmax; i++) {
        for(PS::S64 j = 0; j < nmax; j++) {
            PS::F64 pinv = ((ptcl[i][j] != 0) ? (1. / ptcl[i][j]) : 0.);
            dens[i][j] *= pinv;
            temp[i][j] *= pinv;
            entr[i][j] *= pinv;
        }
    }

    FILE * fp = fopen(ofile, "w");
    for(PS::S64 i = 0; i < nmax; i++) {
        for(PS::S64 j = 0; j < nmax; j++) {
            if(ptcl[i][j] == 0) {
                continue;
            }
            PS::F64 px = xmin[0] + dx * (PS::F64)i;
            PS::F64 py = xmin[1] + dx * (PS::F64)j;
            fprintf(fp, "%+e %+e %+e %+e %+e %8lld\n",
                    px, py,
                    dens[i][j], temp[i][j], entr[i][j],
                    ptcl[i][j]);
        }
    }
    fclose(fp);
}

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);
    PS::ParticleSystem<BlackHoleNeutronStarAnalysis> bhns;
    bhns.initialize();
    bhns.createParticle(0);
    bhns.setNumberOfParticleLocal(0);

    char idir[1024], otype[1024];
    PS::S32 fflag, nfile;
    PS::S64 ibgn, iend;
    PS::F64vec xmin;
    PS::F64    wdth;
    PS::S64    nnxx;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", otype);
    fscanf(fp, "%d %d", &fflag, &nfile);
    fscanf(fp, "%lld%lld", &ibgn, &iend);
    fscanf(fp, "%lf%lf%lf", &xmin[0], &xmin[1], &xmin[2]);
    fscanf(fp, "%lf%lld", &wdth, &nnxx);
    fclose(fp);

    for(PS::S64 itime = ibgn; itime <= iend; itime++) {
        if(fflag == 0) {
            char sfile[1024], bfile[1024];
            sprintf(sfile, "%s/sph_t%04d.dat",  idir, itime);
            sprintf(bfile, "%s/bhns_t%04d.dat", idir, itime);
            fp = fopen(sfile, "r");
            assert(fp);
            sph.readParticleAscii(sfile);
            fclose(fp);
            fp = fopen(bfile, "r");
            assert(fp);
            bhns.readParticleAscii(bfile);
            fclose(fp);
        } else {
            PS::ParticleSystem<SPHAnalysis> tmpsph;
            tmpsph.initialize();
            tmpsph.createParticle(0);
            char sfile[1024];
            PS::F64 ntot = 0;
            for(PS::S64 ifile = 0; ifile < nfile; ifile++) {
                sprintf(sfile, "%s/sph_t%04d_p%06d_i%06d.dat",  idir, itime, nfile, ifile);
                fp = fopen(sfile, "r");
                assert(fp);
                tmpsph.readParticleAscii(sfile);
                fclose(fp);
                PS::S64 ntmp = tmpsph.getNumberOfParticleLocal();
                sph.setNumberOfParticleLocal(ntot+ntmp);
                for(PS::S64 i = 0; i < ntmp; i++) {
                    sph[ntot+i] = tmpsph[i];
                }
                ntot += ntmp;
            }
            char bfile[1024];
            sprintf(bfile, "%s/bhns_t%04d.dat", idir, itime);
            fp = fopen(bfile, "r");
            assert(fp);
            bhns.readParticleAscii(bfile);
            fclose(fp);
        }

        char ofile[1024];
        sprintf(ofile, "%s_t%04d.dat", otype, itime);
        projectOnOrbitalPlane(ofile, sph, bhns, xmin, wdth, nnxx);
    }

    MPI_Finalize();

    return 0;
}
