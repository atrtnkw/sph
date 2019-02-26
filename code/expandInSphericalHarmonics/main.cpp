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
#ifdef USE_INTRINSICS
#include "vector_x86.hpp"
#endif
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

};

namespace SphericalHarmonics {
    PS::F64 CoefficientY00 = 0.5 * sqrt(1.0 / M_PI);
    PS::F64 CoefficientY10 = 0.5 * sqrt(3.0 / M_PI);
    PS::F64 CoefficientY11 = 0.5 * sqrt(1.5 / M_PI);

    bool    FirstTime         = true;
    PS::F64 MinimumRadius     = 1e+08;
    PS::F64 MaximumRadius     = 1e+11;
    PS::S64 NumberOfBinPerDex = 10;

    PS::S64 NumberOfBin = 0;
    const PS::S64 MaximumNumberOfBin = 64;
    static PS::F64 m00re[MaximumNumberOfBin];
    static PS::F64 m10re[MaximumNumberOfBin];
    static PS::F64 m11re[MaximumNumberOfBin];
    static PS::F64 m11im[MaximumNumberOfBin];
    static PS::S64 nptcl[MaximumNumberOfBin];
    static PS::F64 m00re_g[MaximumNumberOfBin];
    static PS::F64 m10re_g[MaximumNumberOfBin];
    static PS::F64 m11re_g[MaximumNumberOfBin];
    static PS::F64 m11im_g[MaximumNumberOfBin];
    static PS::S64 nptcl_g[MaximumNumberOfBin];

    void initialize() {
        FirstTime = false;
        PS::F64 temp = log10(MaximumRadius / MinimumRadius) * NumberOfBinPerDex;
        NumberOfBin = (temp - (PS::S64)temp == 0.) ? (PS::S64)temp : (PS::S64)(temp + 1.);
        assert(NumberOfBin <= MaximumNumberOfBin);
        for(PS::S64 i = 0; i < NumberOfBin; i++) {
            m00re[i] = 0.;
            m10re[i] = 0.;
            m11re[i] = 0.;
            m11im[i] = 0.;
            nptcl[i] = 0;
            m00re_g[i] = 0.;
            m10re_g[i] = 0.;
            m11re_g[i] = 0.;
            m11im_g[i] = 0.;
            nptcl_g[i] = 0.;
        }
    }

    void reinitialize() {
        for(PS::S64 i = 0; i < NumberOfBin; i++) {
            m00re[i] = 0.;
            m10re[i] = 0.;
            m11re[i] = 0.;
            m11im[i] = 0.;
            nptcl[i] = 0;
            m00re_g[i] = 0.;
            m10re_g[i] = 0.;
            m11re_g[i] = 0.;
            m11im_g[i] = 0.;
            nptcl_g[i] = 0.;
        }
    }

    template <class Tsph>
    inline void addParticle(Tsph & isph) {
        PS::F64 r2 = isph.pos * isph.pos;
        if(r2 >= MaximumRadius * MaximumRadius) {
            return;
        }
        PS::F64 r1 = sqrt(r2);
        PS::S64 ibin = PS::S64(log10(r1 / MinimumRadius) * NumberOfBinPerDex);
        assert(ibin < NumberOfBin);

        PS::F64 cr = sqrt(isph.pos[0] * isph.pos[0] + isph.pos[1] * isph.pos[1]);
        PS::F64 prinv = 1. / r1;
        PS::F64 crinv = 1. / cr;
        PS::F64 costheta = isph.pos[2] * prinv;
        PS::F64 sintheta = cr          * prinv;
        PS::F64 cosphi   = isph.pos[0] * crinv;
        PS::F64 sinphi   = isph.pos[1] * crinv;

        m00re[ibin] += isph.mass * CoefficientY00;
        m10re[ibin] += isph.mass * CoefficientY10 * costheta;
        m11re[ibin] += isph.mass * CoefficientY11 * sintheta * cosphi;
        m11im[ibin] += isph.mass * CoefficientY11 * sintheta * sinphi;
        nptcl[ibin] += 1;
    }

    void reduce() {
        PS::S64 ierr = 0;
        ierr = MPI_Allreduce(m00re, m00re_g, NumberOfBin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        ierr = MPI_Allreduce(m10re, m10re_g, NumberOfBin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        ierr = MPI_Allreduce(m11re, m11re_g, NumberOfBin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        ierr = MPI_Allreduce(m11im, m11im_g, NumberOfBin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        ierr = MPI_Allreduce(nptcl, nptcl_g, NumberOfBin, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    }

    void output(FILE * fp) {
        for(PS::S64 i = 0; i < NumberOfBin; i++) {
            PS::F64 r = (i == 0) ? (0.5 * MinimumRadius)
                : (MinimumRadius * pow(10., (i - 0.5) / NumberOfBinPerDex));
            fprintf(fp, "%+e %+e %+e %+e %+e %8d\n",
                    r, m00re_g[i], m10re_g[i], m11re_g[i], m11im_g[i], nptcl_g[i]);
        }
    }

};

template <class Tsph>
void calcDensityCenter(Tsph  & sph,
                       PS::F64vec & xdens_g,
                       PS::F64vec & vdens_g,
                       PS::F64vec & adens_g) {
    PS::F64    tdens = 0.;
    PS::F64vec xdens = 0.;
    PS::F64vec vdens = 0.;
    PS::F64vec adens = 0.;

    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        if(sph[i].istar == 1) {
            continue;
        }
        tdens += sph[i].dens;
        xdens += sph[i].dens * sph[i].pos;
        vdens += sph[i].dens * sph[i].vel;
        adens += sph[i].dens * sph[i].acc;
    }

    PS::F64    tdens_g = PS::Comm::getSum(tdens);
    xdens_g = PS::Comm::getSum(xdens);
    vdens_g = PS::Comm::getSum(vdens);
    adens_g = PS::Comm::getSum(adens);

    xdens_g = xdens_g / tdens_g;
    vdens_g = vdens_g / tdens_g;
    adens_g = adens_g / tdens_g;

    return;
}

template <class Tsph>
void expandInSphericalHarmonics(char * ofile,
                                Tsph  & sph) {
    PS::F64vec xdens = 0.;
    PS::F64vec vdens = 0.;
    PS::F64vec adens = 0.;
    calcDensityCenter(sph, xdens, vdens, adens);
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].pos -= xdens;
        sph[i].vel -= vdens;
        sph[i].acc -= adens;
    }

    if(SphericalHarmonics::FirstTime) {
        SphericalHarmonics::initialize();
    } else {
        SphericalHarmonics::reinitialize();
    }

    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        SphericalHarmonics::addParticle(sph[i]);
    }

    SphericalHarmonics::reduce();

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        SphericalHarmonics::output(fp);
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
    PS::S64 ibgn, iend, dsnp;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld%lld%lld", &ibgn, &iend, &dsnp);
    fclose(fp);

    for(PS::S64 itime = ibgn; itime <= iend; itime += dsnp) {        
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
        sprintf(ofile, "%s/harmonics_t%04d.dat", odir, itime);
        expandInSphericalHarmonics(ofile, sph);
    }

    PS::Finalize();

    return 0;
}
