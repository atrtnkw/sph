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
    const PS::S64 NumberOfMesh    = 200;
    const PS::S64 NumberOfMesh2   = NumberOfMesh * NumberOfMesh;
    const PS::F64 HalfLengthOfBox = 5e8;
    const PS::F64 WidthOfMesh     = (2. * HalfLengthOfBox) / ((PS::F64)NumberOfMesh);
    const PS::F64 WidthOfMeshInv  = 1. / WidthOfMesh;
    static PS::F64 dens_l[NumberOfMesh][NumberOfMesh];
    static PS::F64 dens_g[NumberOfMesh][NumberOfMesh];
    static PS::F64 velz_l[NumberOfMesh][NumberOfMesh];
    static PS::F64 velz_g[NumberOfMesh][NumberOfMesh];
    static PS::F64 velp_l[NumberOfMesh][NumberOfMesh];
    static PS::F64 velp_g[NumberOfMesh][NumberOfMesh];
    static PS::F64 velv_l[NumberOfMesh][NumberOfMesh];
    static PS::F64 velv_g[NumberOfMesh][NumberOfMesh];
    
    const PS::F64vec BasePosition0(+5.854266e+08, -1.250326e+09, 0.);
    const PS::F64vec BasePosition1(+1.594931e+08, -1.375774e+09, 0.);
    PS::F64vec BasisVector;
    PS::F64vec MiddlePoint;
    PS::F64vec StartPoint;

    void initializeAssignToMesh() {
        PS::F64vec dr = BasePosition0 - BasePosition1;
        BasisVector = (1. / sqrt(dr * dr)) * dr;
        MiddlePoint = 0.5 * (BasePosition0 + BasePosition1);
        StartPoint  = MiddlePoint - HalfLengthOfBox * BasisVector;
    }

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

    initializeAssignToMesh();

    for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
        for(PS::S64 is = 0; is < NumberOfMesh; is++) {
            dens_l[iz][is] = 0.;
            velz_l[iz][is] = 0.;
            velp_l[iz][is] = 0.;
            velv_l[iz][is] = 0.;
        }
    }

    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        PS::F64vec pvec2d(sph[ip].pos[0], sph[ip].pos[1], 0.);
        PS::F64    ksr2    = sph[ip].ksr * sph[ip].ksr;
        PS::F64vec dstart  = pvec2d - StartPoint;
        PS::F64    dstart2 = dstart * dstart;
        PS::F64vec dend    = pvec2d - (StartPoint + 2. * HalfLengthOfBox * BasisVector);
        PS::F64    dend2   = dend * dend;
        PS::F64    project = dstart * BasisVector;
        PS::F64    foot2   = dstart2 - project * project;
        if(foot2 > ksr2
           || (project < 0.                   && dstart2 > ksr2)
           || (project > 2. * HalfLengthOfBox && dend2   > ksr2)) {
            continue;
        }
        PS::F64vec vvec2d(sph[ip].vel[0], sph[ip].vel[1], 0.);

        for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
            PS::F64 posz = getPosition(iz);
            for(PS::S64 is = 0; is < NumberOfMesh; is++) {
                PS::F64 length = getPosition(is);
                PS::F64vec pos = MiddlePoint + length * BasisVector;
                pos[2] = posz;

                PS::F64vec dr = pos - sph[ip].pos;
                PS::F64 dr2   = dr * dr;
                PS::F64 dr1   = sqrt(dr2);
                PS::F64 hinv1 = 1. / sph[ip].ksr;
                PS::F64 hinv3 = ND::calcVolumeInverse(hinv1);
                PS::F64 dq1   = dr1 * hinv1;
                PS::F64 kw0   = hinv3 * SK::kernel0th(dq1);
                dens_l[iz][is] += sph[ip].mass * kw0;
                PS::F64 dinv = 1. / sph[ip].dens;
                PS::F64 ksph = sph[ip].mass * dinv * kw0;
                velz_l[iz][is] += ksph * sph[ip].vel[2];
                PS::F64 velp = vvec2d * BasisVector;
                velp_l[iz][is] += ksph * velp;
                PS::F64vec vtemp = vvec2d - velp * BasisVector;
                PS::F64    velv  = sqrt(vtemp * vtemp);
                velv_l[iz][is] += ksph * velv;
            }
        }
    }

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(dens_l, dens_g, NumberOfMesh2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(velz_l, velz_g, NumberOfMesh2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(velp_l, velp_g, NumberOfMesh2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(velv_l, velv_g, NumberOfMesh2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        PS::S64 isstd   = getIndex(0.);
        PS::S64 izstd   = getIndex(0.);
        PS::F64 velpstd = velp_g[izstd][isstd];
        PS::F64 velvstd = velv_g[izstd][isstd];
        for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
            PS::F64 posz = getPosition(iz);
            for(PS::S64 is = 0; is < NumberOfMesh; is++) {
                PS::F64 length = getPosition(is);
                PS::F64vec pos = MiddlePoint + length * BasisVector;
                pos[2] = posz;
                PS::F64 dens = (dens_g[iz][is] > 1e-6) ? dens_g[iz][is] : 1e-6;
                fprintf(fp, "%+e %+e %+e %+e %+e %+e %+e %+e\n",
                        length, pos[0], pos[1], pos[2], dens,
                        velz_g[iz][is], velp_g[iz][is] - velpstd, velv_g[iz][is] - velvstd);
            }
        }
        fclose(fp);
        
        {
            FILE * fp = fopen("t0102.dat.200", "w");
            PS::F64vec relvec  = BasePosition0 - MiddlePoint;
            PS::F64    relvec1 = sqrt(relvec * relvec);
            PS::S64    relis   = getIndex(relvec1);
            for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
                PS::F64 posz = getPosition(iz);
                PS::F64 length = getPosition(relis);
                PS::F64vec pos = MiddlePoint + length * BasisVector;
                pos[2] = posz;
                PS::F64 dens = (dens_g[iz][relis] > 1e-6) ? dens_g[iz][relis] : 1e-6;
                //fprintf(fp, "%+e %+e %+e %+e %+e %+e\n", length, pos[0], pos[1], pos[2], dens, velz_g[iz][relis]);
                fprintf(fp, "%+e %+e %+e\n", pos[2], dens, velz_g[iz][relis]);
            }
            fclose(fp);
        }

        {
            FILE * fp = fopen("t0103.dat.200", "w");
            PS::F64vec relvec  = BasePosition1 - MiddlePoint;
            PS::F64    relvec1 = - sqrt(relvec * relvec);
            PS::S64    relis   = getIndex(relvec1);
            for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
                PS::F64 posz = getPosition(iz);
                PS::F64 length = getPosition(relis);
                PS::F64vec pos = MiddlePoint + length * BasisVector;
                pos[2] = posz;
                PS::F64 dens = (dens_g[iz][relis] > 1e-6) ? dens_g[iz][relis] : 1e-6;
                //fprintf(fp, "%+e %+e %+e %+e %+e %+e\n", length, pos[0], pos[1], pos[2], dens, velz_g[iz][relis]);
                fprintf(fp, "%+e %+e %+e\n", pos[2], dens, velz_g[iz][relis]);
            }
            fclose(fp);        
        }
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
