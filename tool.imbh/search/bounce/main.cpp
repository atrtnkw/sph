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
        fprintf(fp, " %+e", this->pot3);                                       // 45
        fprintf(fp, " %+e %+e %+e", this->tempmax[0], this->tempmax[1], this->tempmax[2]);
                                                                               // 46 -- 48
        fprintf(fp, " %+e", this->entr);                                       // 49
        fprintf(fp, "\n");
    }
};

template <class Tsph>
void search1Dpoint(Tsph & sph,
                   const PS::S64 nx,
                   const PS::S64 ny,
                   const PS::F64 bx,
                   const PS::F64 by,
                   const PS::F64 wx,
                   PS::S64 * nptcl,
                   PS::S64 * pdivv,
                   PS::F64 * dnmax,
                   bool maxornot,
                   bool rightornot,
                   PS::F64 dnlimit,
                   PS::F64 pxlimit,
                   char * filename,
                   PS::F64 * _dnglb,
                   PS::F64 * _pxglb) {
    PS::F64 dnloc = -1.;
    PS::S64 idloc = -1;
    for(PS::S64 i = 0; i < nx * ny; i++) {
        PS::F64 fracpdivv = (PS::F64)pdivv[i] / (PS::F64)nptcl[i];
        if(fracpdivv < 0.3 || 0.6 < fracpdivv) {
            continue;
        }
        if(dnmax[i] > dnlimit && (!maxornot)) {
            continue;
        }
        if(bx + (PS::F64)(i % nx) * wx - pxlimit > 0. && (!maxornot) && (!rightornot)) {
            continue;
        }
        if(dnloc < dnmax[i]) {
            dnloc = dnmax[i];
            idloc = i;
        }
    }

    PS::F64 dnglb = -1.;
    PS::S32 rkglb = -1;
    PS::Comm::getMaxValue(dnloc, PS::Comm::getRank(), dnglb, rkglb);
    PS::S64 idglb = ((PS::Comm::getRank() == rkglb) ? idloc : -1);
    PS::F64 pxglb = bx + (PS::F64)(idloc % nx) * wx;
    PS::Comm::broadcast(&pxglb, 1, rkglb);

    if(PS::Comm::getRank() == rkglb) {
        char ofile[1024];
        sprintf(ofile, filename);
        FILE * fp = fopen(ofile, "w");
        assert(fp);
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            PS::F64 dx = sph[i].pos[0] - bx;
            PS::F64 dy = sph[i].pos[1] - by;
            PS::S64 id = 0;
            id  = (PS::S64)(dy / wx) * nx;
            id += (PS::S64)(dx / wx) % nx;
            assert(id < nx * ny);
            if(id == idglb) {
                sph[i].writeAscii(fp);
            }
        }
        fclose(fp);
    }

    *_dnglb = dnglb;
    *_pxglb = pxglb;
}
    
int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);

    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.setDomain(PS::Comm::getNumberOfProc(), 1);
    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);

    char itype[1024];
    PS::S32 fflag, nfile;
    PS::S64 ibgn, iend;
    PS::F64vec xmin;
    PS::F64    xmax;
    PS::F64    wdth;
    PS::S64    nnxx;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", itype);
    fscanf(fp, "%d %d", &fflag, &nfile);
    fscanf(fp, "%lf%lf%lf", &xmin[0], &xmin[1], &xmax);
    fscanf(fp, "%lf%lld", &wdth, &nnxx);
    fclose(fp);

    if(PS::Comm::getRank() == 0) {
        fprintf(stderr, "itype: %s\n", itype);
        fprintf(stderr, "fflag: %d nfile: %d\n", fflag, nfile);
        fprintf(stderr, "xmin[0]: %+e xmin[1]: %+e\n", xmin[0], xmin[1]);
        fprintf(stderr, "xmax: %+e\n", xmax);
        fprintf(stderr, "width: %+e nnxx: %lld\n", wdth, nnxx);
    }

    // input
    if(fflag == 0) {
        char sfile[1024];
        sprintf(sfile, "%s.dat",  itype);
        fp = fopen(sfile, "r");
        assert(fp);
        fprintf(stderr, "Reading file...\n");
        sph.readParticleAscii(sfile);
        fclose(fp);
    } else {
        PS::S64 ntot = 0;
        PS::S64 nrank = PS::Comm::getNumberOfProc();
        PS::S64 irank = PS::Comm::getRank();
        PS::S64 ihead = nfile *  irank      / nrank;
        PS::S64 itail = nfile * (irank + 1) / nrank;
        for(PS::S64 ifile = ihead; ifile < itail; ifile++) {
            char sfile[1024];
            sprintf(sfile, "%s_p%06d_i%06d.dat",  itype, nfile, ifile);
            fp = fopen(sfile, "r");
            assert(fp);
            PS::S64 ntmp = 0;
            for(PS::S32 c; (c = getc(fp)) != EOF; ntmp += ('\n' == c ? 1 : 0)) {
                ;
            }
            fclose(fp);
            sph.setNumberOfParticleLocal(ntot+ntmp);
            fp = fopen(sfile, "r");
            for(PS::S64 i = 0; i < ntmp; i++) {
                sph[ntot+i].readAscii(fp);
            }
            fclose(fp);
            ntot += ntmp;
        }
    }

    // Data distribution
    {
        PS::S64 nrank = PS::Comm::getNumberOfProc();
        PS::F64 dx    = wdth / (PS::F64)nrank;
        for(PS::S64 i = 0; i < nrank; i++) {
            PS::F64ort pos;
            pos.low_[0]  = xmin[0] + dx * i;
            pos.low_[1]  = - xmax;
            pos.low_[2]  = - xmax;
            pos.high_[0] = xmin[0] + dx * (i + 1);
            pos.high_[1] = + xmax;
            pos.high_[2] = + xmax;
            dinfo.setPosDomain(i, pos);
        }
        sph.exchangeParticle(dinfo);
    }

    // Analyse divergence v
    {
        const PS::F64 wx = wdth / (PS::F64)nnxx;
        const PS::F64 bx = xmin[0] + (wdth / (PS::F64)PS::Comm::getNumberOfProc())
            * PS::Comm::getRank();
        const PS::F64 by = xmin[1];
        const PS::S64 nx = nnxx / PS::Comm::getNumberOfProc();
        const PS::S64 ny = nnxx;
        PS::S64 *nptcl = (PS::S64 *)malloc(sizeof(PS::S64) * nx * ny);
        assert(nptcl);
        PS::S64 *pdivv = (PS::S64 *)malloc(sizeof(PS::S64) * nx * ny);
        assert(pdivv);
        PS::F64 *dnmax = (PS::F64 *)malloc(sizeof(PS::F64) * nx * ny);
        assert(dnmax);
        for(PS::S64 i = 0; i < nx * ny; i++) {
            nptcl[i] = 0;
            pdivv[i] = 0;
            dnmax[i] = 0.;
        }        
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            PS::F64 dx = sph[i].pos[0] - bx;
            PS::F64 dy = sph[i].pos[1] - by;
            PS::S64 id = 0;
            id  = (PS::S64)(dy / wx) * nx;
            id += (PS::S64)(dx / wx) % nx;
            assert(id < nx * ny);
            nptcl[id]++;
            if(sph[i].divv > 0.) {
                pdivv[id]++;
            }
            if(dnmax[id] < sph[i].dens) {
                dnmax[id] = sph[i].dens;
            }
        }

        PS::F64 dncen, pxcen;
        search1Dpoint(sph, nx, ny, bx, by, wx, nptcl, pdivv, dnmax,
                      true, false, 0., 0., "fuga_max.dat", &dncen, &pxcen);
        PS::F64 dntmp, pxtmp;
        search1Dpoint(sph, nx, ny, bx, by, wx, nptcl, pdivv, dnmax,
                      false, true,  0.5*dncen, pxcen, "fuga_rght.dat", &dntmp, &pxtmp);
        search1Dpoint(sph, nx, ny, bx, by, wx, nptcl, pdivv, dnmax,
                      false, false, 0.5*dncen, pxcen, "fuga_left.dat", &dntmp, &pxtmp);

        free(nptcl);
        free(pdivv);
        free(dnmax);
    }

    PS::Finalize();

    return 0;
}
