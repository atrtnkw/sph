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

bool searchIdentity(PS::S64 idi,
                    PS::S64 nlist,
                    PS::S64 * ilist) {
    bool flag         = false;
    PS::S64 imin      = 0;
    PS::S64 imax      = nlist - 1;
    PS::S64 ncount    = 0;
    PS::S64 ncountmax = (1L << 10L);
    do {
        PS::S64 imid = (imin + imax) / 2;
        if(idi == ilist[imin] || idi == ilist[imax]) {
            flag = true;
            break;
        } else if(idi < ilist[imid]) {
            if(imax == imid) break;
            imax = imid;
        } else {
            if(imin == imid) break;
            imin = imid;
        }
        ncount++;
        assert(ncount < ncountmax);
    } while(imax != imin);

    return flag;
}

template <class Tsph>
void searchParticle(char * ofile,
                    Tsph & sph,
                    PS::S64 nlist,
                    PS::S64 * ilist) {
    FILE * fp = fopen(ofile, "w");
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        if(searchIdentity(sph[i].id, nlist, ilist)) {
            sph[i].writeAscii(fp);
        }
    }
    fclose(fp);
}

int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);

    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);

    char idir[1024], odir[1024], lfile[1024];
    PS::S64 ibgn, iend;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld%lld", &ibgn, &iend);
    fscanf(fp, "%s", lfile);
    fclose(fp);

    PS::S64 nlist;
    PS::S64 * ilist;
    fp = fopen(lfile, "r");
    fscanf(fp, "%lld", &nlist);
    ilist = (PS::S64 *)malloc(sizeof(PS::S64) * nlist);
    assert(ilist);
    for(PS::S64 i = 0; i < nlist; i++) {
        fscanf(fp, "%lld", &ilist[i]);
    }
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
            continue;
        }
        fclose(fp);

        char sfile[1024];
        sprintf(sfile, "%s/t%02d/sph_t%04d", idir, tdir, itime);
        sph.readParticleAscii(sfile, "%s_p%06d_i%06d.dat");

        char ofile[1024];
        sprintf(ofile, "%s/sph_t%04d_p%06d_i%06d.dat",
                odir, itime, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
        searchParticle(ofile, sph, nlist, ilist);

    }

    free(ilist);

    PS::Finalize();

    return 0;
}
