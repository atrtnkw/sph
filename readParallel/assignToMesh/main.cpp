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
#define SurvivingWhiteDwarf
#define Data2D
#define LowResolution
#ifdef LowResolution
    const PS::S64 NumberOfMesh    = 128;
#else
    const PS::S64 NumberOfMesh    = 256;
#endif
    const PS::S64 NumberOfMesh2   = NumberOfMesh * NumberOfMesh;
    const PS::S64 NumberOfMesh3   = NumberOfMesh * NumberOfMesh2;
#ifdef SurvivingWhiteDwarf
    const PS::F64 HalfLengthOfBox = 2.e9;
#else
    const PS::F64 HalfLengthOfBox = 3.e11;
#endif
    const PS::F64 WidthOfMesh     = (2. * HalfLengthOfBox) / ((PS::F64)NumberOfMesh);
    const PS::F64 WidthOfMeshInv  = 1. / WidthOfMesh;
    static PS::F32 dens_l[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 dens_g[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 uene_l[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 uene_g[NumberOfMesh][NumberOfMesh][NumberOfMesh];
#ifdef LowResolution
    static PS::F32 temp_l[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 temp_g[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 entr_l[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 entr_g[NumberOfMesh][NumberOfMesh][NumberOfMesh];
#endif
    static PS::F32 cmps_l[NR::NumberOfNucleon][NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 cmps_g[NR::NumberOfNucleon][NumberOfMesh][NumberOfMesh][NumberOfMesh];

    PS::S64 getIndex(PS::F64 x) {
        PS::F64 fx = (x + HalfLengthOfBox) * WidthOfMeshInv - 0.5;
        PS::S64 ix = (PS::S64)(ceil(fx));
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
void assignToMesh(char * ofile,
                  Tsph & sph) {
    using namespace AssignToMesh;

    PS::F64vec pcenter(0.);
    PS::F64vec vcenter(0.);
    obtainDensityCenter(sph, pcenter, vcenter);

    for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
        for(PS::S64 iy = 0; iy < NumberOfMesh; iy++) {
            for(PS::S64 ix = 0; ix < NumberOfMesh; ix++) {
                dens_l[iz][iy][ix] = 0.;
                uene_l[iz][iy][ix] = 0.;
#ifdef LowResolution
                temp_l[iz][iy][ix] = 0.;
                entr_l[iz][iy][ix] = 0.;
#endif
                for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                    cmps_l[k][iz][iy][ix] = 0.;
                }
            }
        }
    }

    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        PS::F64vec pptcl = sph[ip].pos - pcenter;
        PS::F64vec vptcl = sph[ip].vel - vcenter;
        PS::F64    eptcl = 0.5 * (vptcl * vptcl) + sph[ip].pot;
#ifdef SurvivingWhiteDwarf
        if(eptcl >= 0.) {
            continue;
        }
#else
        if(eptcl < 0.) {
            continue;
        }
#endif

        PS::S64 ixmin = getIndex(pptcl[0] - sph[ip].ksr);
        PS::S64 ixmax = getIndex(pptcl[0] + sph[ip].ksr);
        PS::S64 iymin = getIndex(pptcl[1] - sph[ip].ksr);
        PS::S64 iymax = getIndex(pptcl[1] + sph[ip].ksr);
        PS::S64 izmin = getIndex(pptcl[2] - sph[ip].ksr);
        PS::S64 izmax = getIndex(pptcl[2] + sph[ip].ksr);
        assert(0 <= ixmin && ixmin < NumberOfMesh);
        assert(0 <= ixmax && ixmax < NumberOfMesh);
        assert(0 <= iymin && iymin < NumberOfMesh);
        assert(0 <= iymax && iymax < NumberOfMesh);
        assert(0 <= izmin && izmin < NumberOfMesh);
        assert(0 <= izmax && izmax < NumberOfMesh);
//        for(PS::S64 iz = izmin; iz < izmax; iz++) {
//            for(PS::S64 iy = iymin; iy < iymax; iy++) {
//                for(PS::S64 ix = ixmin; ix < ixmax; ix++) {
        for(PS::S64 iz = izmin; iz <= izmax; iz++) {
            for(PS::S64 iy = iymin; iy <= iymax; iy++) {
                for(PS::S64 ix = ixmin; ix <= ixmax; ix++) {
                    PS::F64vec pos(getPosition(ix),
                                   getPosition(iy),
                                   getPosition(iz));
                    PS::F64vec dr = pos - pptcl;
                    PS::F64 dr2   = dr * dr;
                    PS::F64 dr1   = sqrt(dr2);
                    PS::F64 hinv1 = 1. / sph[ip].ksr;
                    PS::F64 hinv3 = ND::calcVolumeInverse(hinv1);
                    PS::F64 dq1   = dr1 * hinv1;
                    PS::F64 kw0   = hinv3 * SK::kernel0th(dq1);
                    dens_l[iz][iy][ix] += sph[ip].mass * kw0;

                    PS::F64 dinv = 1. / sph[ip].dens;
                    PS::F64 ksph = sph[ip].mass * dinv * kw0;
                    uene_l[iz][iy][ix] += ksph * sph[ip].uene;
#ifdef LowResolution
                    temp_l[iz][iy][ix] += ksph * sph[ip].temp;
                    entr_l[iz][iy][ix] += ksph * sph[ip].entr;
#endif
                    for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                        cmps_l[k][iz][iy][ix] += ksph * sph[ip].cmps[k];
                    }

                }            
            }            
        }
    }

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(dens_l, dens_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(uene_l, uene_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#ifdef LowResolution
    ierr = MPI_Allreduce(temp_l, temp_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(entr_l, entr_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
    for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
        ierr = MPI_Allreduce(cmps_l[k], cmps_g[k], NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }

#ifdef Data2D
    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        for(PS::S64 iy = 0; iy < NumberOfMesh; iy++) {
            for(PS::S64 ix = 0; ix < NumberOfMesh; ix++) {
                PS::F64 mx = getPosition(ix);
                PS::F64 my = getPosition(iy);
                fprintf(fp, "%+e %+e %+e %+e", mx, my,
                        dens_g[NumberOfMesh/2][iy][ix],
                        uene_g[NumberOfMesh/2][iy][ix]);
#ifdef LowResolution
                fprintf(fp, " %+e", temp_g[NumberOfMesh/2][iy][ix]);
                fprintf(fp, " %+e", entr_g[NumberOfMesh/2][iy][ix]);
#endif
                PS::F32 norm = 0.;
                for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                    norm += cmps_g[k][NumberOfMesh/2][iy][ix];
                }
                PS::F32 ninv = (norm != 0.) ? (1. / norm) : 0.;
                for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                    fprintf(fp, " %+.2e", cmps_g[k][NumberOfMesh/2][iy][ix] * ninv);
                }
                fprintf(fp, "\n");
            }
        }
        fclose(fp);
    }
#else
    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
            for(PS::S64 iy = 0; iy < NumberOfMesh; iy++) {
                for(PS::S64 ix = 0; ix < NumberOfMesh; ix++) {
                    PS::F64 mx = getPosition(ix);
                    PS::F64 my = getPosition(iy);
                    PS::F64 mz = getPosition(iz);
#if 0
                    fprintf(fp, "%+e %+e %+e %+e", mx, my, mz, dens_g[iz][iy][ix]);
#else
                    PS::F32 norm = 0.;
                    for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                        norm += cmps_g[k][NumberOfMesh/2][iy][ix];
                    }
                    PS::F32 ninv = (norm != 0.) ? (1. / norm) : 0.;
                    fprintf(fp, " %+e %+e %+e", mx, my, mz);
                    fprintf(fp, " %+e %+e %+e %+e", dens_g[iz][iy][ix], uene_g[iz][iy][ix], temp_g[iz][iy][ix], entr_g[iz][iy][ix]);
                    for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                        fprintf(fp, " %+.2e", cmps_g[k][NumberOfMesh/2][iy][ix] * ninv);
                    }
#endif
                    fprintf(fp, "\n");
                }
            }
        }
        fclose(fp);
    }
#endif

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
#ifdef Data2D
        sprintf(ofile, "%s/init2d_t%04d.dat", odir, itime);
#else
        sprintf(ofile, "%s/init3d_t%04d.dat", odir, itime);
#endif
        assignToMesh(ofile, sph);

    }

    PS::Finalize();

    return 0;
}
