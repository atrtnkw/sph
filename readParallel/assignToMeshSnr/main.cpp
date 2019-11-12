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

class SPHDetailElement : public HelmholtzGas {
public:

    static       PS::S64 CompanionCO;
    static const PS::S64 NumberOfElement = 34;
    PS::F64 elem[NumberOfElement];

    SPHDetailElement() {
        this->id    = 0;
        this->istar = 0;
        this->mass  = 0.;
        this->pos   = 0.;
        this->vel   = 0.;
        this->dens  = 0.;
        this->ksr   = 0.;
        this->pot   = 0.;
        for(PS::S64 k = 0; k < NumberOfElement; k++) {
            this->elem[k] = 0.;
        }
    }

    void readAscii(FILE * fp) {
        fscanf(fp, "%lld%lld%lf", &this->id,     &this->istar,  &this->mass);   //   3
        fscanf(fp, "%lf%lf%lf",   &this->pos[0], &this->pos[1], &this->pos[2]); //   6
        fscanf(fp, "%lf%lf%lf",   &this->vel[0], &this->vel[1], &this->vel[2]); //   9
        fscanf(fp, "%lf%lf%lf",   &this->dens,   &this->ksr,    &this->pot);    //  12
        if(CompanionCO == 1) {
            this->elem[9]  = 0.49179301803663630124; // C
            this->elem[11] = 0.49179301803663630124; // O
            this->elem[13] = 0.01342153152708231446; // Ne
            this->elem[15] = 0.00070522296620528707; // Mg
            this->elem[17] = 0.00066876453890332400; // Si
            this->elem[19] = 0.00031136169559311488; // S
            this->elem[29] = 0.00130708319894335710; // Fe
        } else if(CompanionCO == 0) {
            this->elem[5]  = 0.59232608581213036673; // He
            this->elem[9]  = 0.19671720721465452049; // C
            this->elem[10] = 0.00512458476488597461; // N
            this->elem[11] = 0.19671720721465452049; // O
            this->elem[13] = 0.00075386998319660882 + 0.00536861261083292578; // Ne
            this->elem[15] = 0.00070522296620528707; // Mg
            this->elem[17] = 0.00066876453890332400; // Si
            this->elem[19] = 0.00031136169559311488; // S
            this->elem[29] = 0.00130708319894335710; // Fe
        } else {
            assert(NULL);
        }
        if(this->istar == 0) {
            PS::S64 dum;
            fscanf(fp, "%lld", &dum);
            for(PS::S64 k = 0; k < NumberOfElement; k++) {
                fscanf(fp, "%lf", &this->elem[k]);
            }
        }
    }

};
PS::S64 SPHDetailElement::CompanionCO        = -1;

namespace AssignToMesh {
    bool Flag2D = false;
    const PS::S64 MaxNumberOfElement = 10;
    //const PS::S64 NumberOfMesh    = 256;
    const PS::S64 NumberOfMesh    = 64;
    const PS::S64 NumberOfMesh2   = NumberOfMesh * NumberOfMesh;
    const PS::S64 NumberOfMesh3   = NumberOfMesh * NumberOfMesh2;
    const PS::F64 HalfLengthOfBox = 3.e11;
    const PS::F64 WidthOfMesh     = (2. * HalfLengthOfBox) / ((PS::F64)NumberOfMesh);
    const PS::F64 WidthOfMeshInv  = 1. / WidthOfMesh;
    static PS::F32 dens_l[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 dens_g[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 cmps_l[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 cmps_g[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 ni56_l[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 ni56_g[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 etot_l[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 etot_g[NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 elem_l[MaxNumberOfElement][NumberOfMesh][NumberOfMesh][NumberOfMesh];
    static PS::F32 elem_g[MaxNumberOfElement][NumberOfMesh][NumberOfMesh][NumberOfMesh];
    PS::S64 FirstElement = -1;

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


template <class Tsph,
          class Tsphelem>
void assignToMesh(char * ofile,
                  Tsph & sph,
                  Tsphelem & sphelem) {
    using namespace AssignToMesh;

    PS::F64vec pcenter(0.);
    PS::F64vec vcenter(0.);
    obtainDensityCenter(sph, pcenter, vcenter);

    for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
        for(PS::S64 iy = 0; iy < NumberOfMesh; iy++) {
            for(PS::S64 ix = 0; ix < NumberOfMesh; ix++) {
                dens_l[iz][iy][ix] = 0.;
                cmps_l[iz][iy][ix] = 0.;
                ni56_l[iz][iy][ix] = 0.;
                for(PS::S64 k = 0; k < MaxNumberOfElement; k++) {
                    elem_l[k][iz][iy][ix] = 0.;
                }
            }
        }
    }

    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        PS::F64vec pptcl = sph[ip].pos - pcenter;
        PS::F64vec vptcl = sph[ip].vel - vcenter;
        PS::F64    eptcl = 0.5 * (vptcl * vptcl) + sph[ip].pot;
        if(eptcl < 0.) {
            continue;
        }

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
                    cmps_l[iz][iy][ix] += ksph;
                    ni56_l[iz][iy][ix] += ksph * sph[ip].cmps[12];
                }            
            }            
        }
    }

    for(PS::S64 ip = 0; ip < sphelem.getNumberOfParticleLocal(); ip++) {
        PS::F64vec pptcl = sphelem[ip].pos - pcenter;
        PS::F64vec vptcl = sphelem[ip].vel - vcenter;
        PS::F64    eptcl = 0.5 * (vptcl * vptcl) + sphelem[ip].pot;
        if(eptcl < 0.) {
            continue;
        }
        if(sphelem[ip].id % 320 != 0) {
            continue;
        }

        PS::F64 fexp  = pow(320., 1./3.);
        PS::F64 finv3 = 1. / pow(fexp, 3.);
        PS::F64 iksr  = fexp * sphelem[ip].ksr;
        PS::S64 ixmin = getIndex(pptcl[0] - iksr);
        PS::S64 ixmax = getIndex(pptcl[0] + iksr);
        PS::S64 iymin = getIndex(pptcl[1] - iksr);
        PS::S64 iymax = getIndex(pptcl[1] + iksr);
        PS::S64 izmin = getIndex(pptcl[2] - iksr);
        PS::S64 izmax = getIndex(pptcl[2] + iksr);
        assert(0 <= ixmin && ixmin < NumberOfMesh);
        assert(0 <= ixmax && ixmax < NumberOfMesh);
        assert(0 <= iymin && iymin < NumberOfMesh);
        assert(0 <= iymax && iymax < NumberOfMesh);
        assert(0 <= izmin && izmin < NumberOfMesh);
        assert(0 <= izmax && izmax < NumberOfMesh);

        for(PS::S64 iz = izmin; iz <= izmax; iz++) {
            for(PS::S64 iy = iymin; iy <= iymax; iy++) {
                for(PS::S64 ix = ixmin; ix <= ixmax; ix++) {
                    PS::F64vec pos(getPosition(ix),
                                   getPosition(iy),
                                   getPosition(iz));
                    PS::F64vec dr = pos - pptcl;
                    PS::F64 dr2   = dr * dr;
                    PS::F64 dr1   = sqrt(dr2);
                    PS::F64 hinv1 = 1. / iksr;
                    PS::F64 hinv3 = ND::calcVolumeInverse(hinv1);
                    PS::F64 dq1   = dr1 * hinv1;
                    PS::F64 kw0   = hinv3 * SK::kernel0th(dq1);
                    PS::F64 dinv = 1. / (sphelem[ip].dens * finv3);
                    PS::F64 ksph = sphelem[ip].mass * dinv * kw0;
                    etot_l[iz][iy][ix] += ksph;
                    for(PS::S64 k = 0; k < MaxNumberOfElement; k++) {
                        elem_l[k][iz][iy][ix] += ksph * sphelem[ip].elem[k+FirstElement];
                    }
                }
            }
        }
    }

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(dens_l, dens_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(cmps_l, cmps_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(ni56_l, ni56_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(etot_l, etot_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    for(PS::S64 k = 0; k < MaxNumberOfElement; k++) {
        ierr = MPI_Allreduce(elem_l[k], elem_g[k], NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }

    if(PS::Comm::getRank() == 0) {
        if(Flag2D) {
            FILE * fp = fopen(ofile, "w");
            PS::S64 iz = NumberOfMesh / 2;
            for(PS::S64 iy = 0; iy < NumberOfMesh; iy++) {
                for(PS::S64 ix = 0; ix < NumberOfMesh; ix++) {
                    PS::F64 mx   = getPosition(ix);
                    PS::F64 my   = getPosition(iy);
                    fprintf(fp, " %+e %+e", mx, my);
                    fprintf(fp, " %+e", dens_g[iz][iy][ix]);
                    
                    PS::F64 norm =  cmps_g[iz][iy][ix];
                    PS::F64 ninv = (norm != 0.) ? (1. / norm) : 0.;
                    fprintf(fp, " %+e", ni56_g[iz][iy][ix] * ninv);
                    
                    norm = etot_g[iz][iy][ix];
                    ninv = (norm != 0.) ? (1. / norm) : 0.;
                    ninv = (dens_g[iz][iy][ix] != 0.) ? ninv : 0.;
                    for(PS::S64 k = 0; k < MaxNumberOfElement; k++) {
                        fprintf(fp, " %+e", elem_g[k][iz][iy][ix] * ninv);
                    }
                    
                    fprintf(fp, "\n");
                }
            }
            fclose(fp);
        } else {
            FILE * fp = fopen(ofile, "w");
            for(PS::S64 iz = 0; iz < NumberOfMesh; iz++) {
                for(PS::S64 iy = 0; iy < NumberOfMesh; iy++) {
                    for(PS::S64 ix = 0; ix < NumberOfMesh; ix++) {
                        PS::F64 mx   = getPosition(ix);
                        PS::F64 my   = getPosition(iy);
                        PS::F64 mz   = getPosition(iz);
                        fprintf(fp, " %+e %+e %+e", mx, my, mz);
                        fprintf(fp, " %+e", dens_g[iz][iy][ix]);
                        
                        PS::F64 norm =  cmps_g[iz][iy][ix];
                        PS::F64 ninv = (norm != 0.) ? (1. / norm) : 0.;
                        fprintf(fp, " %+e", ni56_g[iz][iy][ix] * ninv);
                        
                        norm = etot_g[iz][iy][ix];
                        ninv = (norm != 0.) ? (1. / norm) : 0.;
                        ninv = (dens_g[iz][iy][ix] != 0.) ? ninv : 0.;
                        for(PS::S64 k = 0; k < MaxNumberOfElement; k++) {
                            fprintf(fp, " %+e", elem_g[k][iz][iy][ix] * ninv);
                        }
                        
                        fprintf(fp, "\n");
                    }
                }
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

    PS::ParticleSystem<SPHDetailElement> sphelem;
    sphelem.initialize();
    sphelem.createParticle(0);
    sphelem.setNumberOfParticleLocal(0);

    char idir[1024], otype[1024], etype[1024];
    PS::S64 itime;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", otype);
    fscanf(fp, "%s", etype);
    fscanf(fp, "%lld", &itime);
    fscanf(fp, "%lld", &SPHDetailElement::CompanionCO);
    fscanf(fp, "%lld", &AssignToMesh::FirstElement);
    fclose(fp);

    assert(AssignToMesh::FirstElement != -1);

    {        
        char tfile[1024];
        FILE *fp = NULL;
        PS::S64 tdir = 0;
        for(PS::S64 iidir = 0; iidir < 100; iidir++) {
            sprintf(tfile, "%s/t%02d/sph_t%04lld_p%06d_i%06d.dat", idir, iidir, itime,
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

        char efile[1024];
        sphelem.readParticleAscii(etype, "%s_p%06d_i%06d.data");

        char ofile[1024];
        sprintf(ofile, "%sE%02d.dat", otype, AssignToMesh::FirstElement);
        assignToMesh(ofile, sph, sphelem);

    }

    PS::Finalize();

    return 0;
}
