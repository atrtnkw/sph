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

    void writeAscii(FILE *fp) const {
        fprintf(fp, "%6d %2d %+e", this->id, this->istar, this->mass);         //  3
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]); //  6
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]); //  9
        fprintf(fp, " %+e %+e %+e", this->acc[0], this->acc[1], this->acc[2]); // 12
        fprintf(fp, " %+e %+e %+e", this->uene, this->alph, this->alphu);      // 15
        fprintf(fp, " %+e %+e %6d", this->dens, this->ksr,  this->np);         // 18
        fprintf(fp, " %+e %+e %+e", this->vsnd, this->pres, this->temp);       // 21
        fprintf(fp, " %+e %+e %+e", this->divv, this->rotv, this->bswt);       // 24
        fprintf(fp, " %+e %+e %+e", this->pot, this->abar, this->zbar);        // 27
        fprintf(fp, " %+e",         this->enuc);                               // 28
        fprintf(fp, " %+e %+e %+e", this->vsmx, this->udot, this->dnuc);       // 31
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) {        // 32 -- 44
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        if(RP::FlagDamping == 2) {
            PS::F64vec omg = RP::RotationalVelocity;
            PS::F64vec vec = omg ^ pos;
            fprintf(fp, " %+e", this->pot - 0.5 * (vec * vec));                // 45
        }
        fprintf(fp, " %+e", this->pot3);                                       // 45
        fprintf(fp, " %+e %+e %+e", this->tempmax[0], this->tempmax[1],
                this->tempmax[2]);                                             // 46 -- 48
        fprintf(fp, " %+e", this->entr);                                       // 49
        fprintf(fp, "\n");
    }

};

template <class Tsph>
void calcPositionAndVelocityOfDensityCenter(Tsph & sph,
                                            PS::F64vec & xcglb,
                                            PS::F64vec & vcglb) {
    PS::F64    dnloc = 0.;
    PS::F64vec xcloc = 0.;
    PS::F64vec vcloc = 0.;
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        dnloc += sph[i].dens;
        xcloc += sph[i].dens * sph[i].pos;
        vcloc += sph[i].dens * sph[i].vel;
    }
    PS::F64 dnglb;
    dnglb = PS::Comm::getSum(dnloc);
    xcglb = PS::Comm::getSum(xcloc);
    vcglb = PS::Comm::getSum(vcloc);
    xcglb /= dnglb;
    vcglb /= dnglb;
}

template <class Tsph>
void findBoundParticle(char * hfile,
                       char * ofile,
                       PS::S64 itime,
                       PS::S64 searchmode,
                       PS::S64 printmode,
                       PS::S64 excludedID,
                       Tsph & sph) {
    PS::F64vec xc, vc;
    calcPositionAndVelocityOfDensityCenter(sph, xc, vc);

    NR::Nucleon mxloc;
    FILE * fp = NULL;
    if(printmode == 0) {
        fp = fopen(ofile, "w");
    }
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        PS::F64vec dx = sph[i].pos - xc;
        PS::F64vec dv = sph[i].vel - vc;
        PS::F64    eb = 0.5 * (dv * dv) + sph[i].pot;
        eb = (searchmode == 0) ? eb : (- eb);
        //if(eb < 0.) {
        if(eb < 0. && !(sph[i].istar == excludedID)) {
            if(printmode == 0) {
                sph[i].writeAscii(fp);
            }
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                mxloc[k] += sph[i].mass * sph[i].cmps[k];
            }
        }
    }
    if(printmode == 0) {
        fclose(fp);
    }
    NR::Nucleon mxglb;
    for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
        mxglb[k]  = PS::Comm::getSum(mxloc[k]);
        mxglb[k] /= CodeUnit::SolarMass;
    }

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(hfile, "a");
        fprintf(fp, "%8d", itime);                                         // 1
        fprintf(fp, " %+e %+e %+e", xc[0], xc[1], xc[2]);                  // 4
        fprintf(fp, " %+e %+e %+e", vc[0], vc[1], vc[2]);                  // 7
        fprintf(fp, " %+.3e %+.3e %+.3e", mxglb[0], mxglb[1],  mxglb[2]);  // 10
        fprintf(fp, " %+.3e %+.3e %+.3e", mxglb[3], mxglb[4],  mxglb[5]);  // 13
        fprintf(fp, " %+.3e %+.3e %+.3e", mxglb[6], mxglb[7],  mxglb[8]);  // 16
        fprintf(fp, " %+.3e %+.3e %+.3e", mxglb[9], mxglb[10], mxglb[11]); // 19
        fprintf(fp, " %+.3e", mxglb[12]);                                  // 20
        fprintf(fp, "\n");
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
    PS::S64 searchmode;
    PS::S64 printmode;
    PS::S64 excludedID;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld%lld", &ibgn, &iend);
    fscanf(fp, "%lld", &searchmode);
    fscanf(fp, "%lld", &printmode);
    fscanf(fp, "%lld", &excludedID);
    fclose(fp);

    char hfile[1024];
    sprintf(hfile, "%s/summary.log", odir);
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
            continue;
        }
        fclose(fp);

        char sfile[1024];
        sprintf(sfile, "%s/t%02d/sph_t%04d", idir, tdir, itime);
        sph.readParticleAscii(sfile, "%s_p%06d_i%06d.dat");

        char ofile[1024];
        //sprintf(ofile, "%s/sph_t%04d_p%06d_i%06d.dat.bound",
        sprintf(ofile, "%s/sph_t%04d_p%06d_i%06d.dat",
                odir, itime, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
        findBoundParticle(hfile, ofile, itime, searchmode, printmode, excludedID, sph);

    }

    PS::Finalize();

    return 0;
}
