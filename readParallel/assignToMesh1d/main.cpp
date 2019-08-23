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
    const  PS::S64 NumberOfMesh   = 256;
    const  PS::F64 MinimumRadius  = 1.e7;
    const  PS::F64 MaximumRadius  = 2.e9;
    const  PS::F64 LogMin         = log(MinimumRadius);
    const  PS::F64 LogMaxMinInv   = 1. / log(MaximumRadius / MinimumRadius);
    const  PS::F64 MaxMinPower    = pow(MaximumRadius / MinimumRadius, 1. / (PS::F64)NumberOfMesh);
    const  PS::S64 NumberOfVortex = 12;
    const  PS::F64 NumberOfVortexInv = 1. / (PS::F64)NumberOfVortex;
    static PS::F64 dens_l[NumberOfMesh];
    static PS::F64 dens_g[NumberOfMesh];
    static PS::F64 uene_l[NumberOfMesh];
    static PS::F64 uene_g[NumberOfMesh];
    static PS::F64 temp_l[NumberOfMesh];
    static PS::F64 temp_g[NumberOfMesh];
    static PS::F64 entr_l[NumberOfMesh];
    static PS::F64 entr_g[NumberOfMesh];
    static PS::F64 cmps_l[NR::NumberOfNucleon][NumberOfMesh];
    static PS::F64 cmps_g[NR::NumberOfNucleon][NumberOfMesh];

    PS::S64 getIndex(PS::F64 r) {
        PS::F64 fr = (log(r) - LogMin) * LogMaxMinInv * (PS::F64)NumberOfMesh;
        PS::S64 ir = (PS::S64)(ceil(fr));
        ir = (ir >= 0)           ? ir : 0;
        ir = (ir < NumberOfMesh) ? ir : (NumberOfMesh - 1);
        return ir;
    }

    PS::F64 getPosition(PS::S64 index) {
        return MinimumRadius * pow(MaxMinPower, (PS::F64)index);
    }

    PS::F64vec getVortexRegularIcosahedron(PS::S64 index) {
        assert(0 <= index < NumberOfVortex);
        PS::F64vec vex = 0.;
        PS::F64 px, py, pz;
        PS::F64 rxy, phi;
        if(index == 0) {
            px =  0.;
            py =  0.;
            pz = +1.;
        } else if (index == 1) {
            px =  0.;
            py =  0.;
            pz = -1.;
        } else if (index <= 6){
            pz  = + sin(atan(0.5));
            rxy = sqrt(1. - pz * pz);
            phi = (index - 2.) * 0.2 * 2. * M_PI;
            px  = rxy * cos(phi);
            py  = rxy * sin(phi);
        } else {
            pz  = - sin(atan(0.5));
            rxy = sqrt(1. - pz * pz);
            phi = (index - 7 + 0.5) * 0.2 * 2. * M_PI;
            px  = rxy * cos(phi);
            py  = rxy * sin(phi);
        }
        vex[0] = px;
        vex[1] = py;
        vex[2] = pz;
        return vex;
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

    for(PS::S64 ir = 0; ir < NumberOfMesh; ir++) {
        dens_l[ir] = 0.;
        uene_l[ir] = 0.;
        temp_l[ir] = 0.;
        entr_l[ir] = 0.;
        for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
            cmps_l[k][ir] = 0.;
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
        PS::F64 rptcl = sqrt(pptcl * pptcl);
        PS::S64 irmin = getIndex(rptcl - sph[ip].ksr);
        PS::S64 irmax = getIndex(rptcl + sph[ip].ksr);
        assert(0 <= irmin && irmin < NumberOfMesh);
        assert(0 <= irmax && irmax < NumberOfMesh);
        for(PS::S64 ir = irmin; ir <= irmax; ir++) {
            for(PS::S64 ivex = 0; ivex < NumberOfVortex; ivex++) {
                PS::F64vec pos = getPosition(ir) * getVortexRegularIcosahedron(ivex);
                PS::F64vec dr = pos - pptcl;
                PS::F64 dr2   = dr * dr;
                PS::F64 dr1   = sqrt(dr2);
                PS::F64 hinv1 = 1. / sph[ip].ksr;
                PS::F64 hinv3 = ND::calcVolumeInverse(hinv1);
                PS::F64 dq1   = dr1 * hinv1;
                PS::F64 kw0   = hinv3 * SK::kernel0th(dq1);
                dens_l[ir]   += sph[ip].mass * kw0 * NumberOfVortexInv;
                PS::F64 dinv  = 1. / sph[ip].dens;
                PS::F64 ksph  = sph[ip].mass * dinv * kw0;
                uene_l[ir] += ksph * sph[ip].uene * NumberOfVortexInv;
                temp_l[ir] += ksph * sph[ip].temp * NumberOfVortexInv;
                entr_l[ir] += ksph * sph[ip].entr * NumberOfVortexInv;
                for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                    cmps_l[k][ir] += ksph * sph[ip].cmps[k] * NumberOfVortexInv;;
                }
            }
        }
    }

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(dens_l, dens_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(uene_l, uene_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(temp_l, temp_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(entr_l, entr_g, NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
        ierr = MPI_Allreduce(cmps_l[k], cmps_g[k], NumberOfMesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        for(PS::S64 ir = 0; ir < NumberOfMesh; ir++) {
            PS::F64 rad = getPosition(ir);

            PS::F64 norm = 0.;
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                norm += cmps_g[k][ir];
            }
            PS::F64 norminv = 1. / norm;
            NR::Nucleon cmps;
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                cmps[k] = cmps_g[k][ir] * norminv;;
            }

            PS::F64 tin = (temp_g[ir] != 0.) ? temp_g[ir] : 1e9;
            PS::F64 pres, vsnd, temp, entr;
            flash_helmholtz_(&dens_g[ir], &uene_g[ir], &tin, cmps.getPointer(),
                             &pres, &vsnd, &temp, &entr);

            fprintf(fp, " %+e", rad);
            fprintf(fp, " %+e", dens_g[ir]);
            fprintf(fp, " %+e", uene_g[ir]);
            fprintf(fp, " %+e", temp_g[ir]);
            fprintf(fp, " %+e", entr_g[ir]);
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                fprintf(fp, " %+.2e", cmps[k]);
            }
            fprintf(fp, " %+e", temp);
            fprintf(fp, " %+e", entr);
            fprintf(fp, "\n");
            fflush(fp);
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
    fclose(fp);

    init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);

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
        sprintf(ofile, "%s/init1d_t%04d.dat", odir, itime);
        assignToMesh(ofile, sph);

    }

    PS::Finalize();

    return 0;
}
