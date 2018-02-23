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

class Plane {
    PS::F64vec nvec;
    PS::F64    cnst;
    PS::F64vec base;
    PS::F64vec axis[2];
    PS::F64    wdth;

    void normalize() {
        PS::F64 vabs = sqrt(this->nvec * this->nvec);
        this->nvec = this->nvec / vabs;
        this->cnst = this->cnst / vabs;
    }

    void defineAxis() {
        if(this->nvec[0] == 0. && this->nvec[1] == 0. && this->nvec[2] == 1.) {
            this->axis[0] = PS::F64vec(1., 0., 0.);
            this->axis[1] = PS::F64vec(0., 1., 0.);
        } else if (this->nvec[0] == 0. && this->nvec[1] == 1. && this->nvec[2] == 0.) {
            this->axis[0] = PS::F64vec(1., 0., 0.);
            this->axis[1] = PS::F64vec(0., 0., 1.);
        } else {
            assert(NULL);
        }
    }

public:

    Plane() {
        this->nvec    = 0.;
        this->cnst    = 0.;
        this->base    = 0.;
        this->axis[0] = 0.;
        this->axis[1] = 0.;
        this->wdth    = 0.;
    }

    PS::F64 readData(FILE *fp) {
        fscanf(fp, "%lf%lf%lf%lf",
               &this->nvec[0], &this->nvec[1], &this->nvec[2], &this->cnst);
        fscanf(fp, "%lf", &this->wdth);
        this->normalize();
        this->base = this->cnst * this->nvec;
        this->defineAxis();
    }

    PS::F64 writeData(FILE *fp) {
        if(PS::Comm::getRank() == 0) {
            fprintf(fp, "nvec: %+e %+e %+e\n",
                    this->nvec[0], this->nvec[1], this->nvec[2]);
            fprintf(fp, "cnst: %+e\n", this->cnst);
            fprintf(fp, "base: %+e %+e %+e\n",
                    this->base[0], this->base[1], this->base[2]);
            fprintf(fp, "axis[0]: %+e %+e %+e\n",
                    this->axis[0][0], this->axis[0][1], this->axis[0][2]);
            fprintf(fp, "axis[1]: %+e %+e %+e\n",
                    this->axis[1][0], this->axis[1][1], this->axis[1][2]);
        }
    }

    inline PS::F64 calcDistance(PS::F64vec pos) {
        return fabs(this->nvec * (pos - this->base));
    }

    inline PS::F64vec calcCoordinateOnPlane(PS::F64vec pos) {
        PS::F64 px = this->axis[0] * (pos - this->base);
        PS::F64 py = this->axis[1] * (pos - this->base);
        return PS::F64vec(px, py, 0.);
    }

    PS::F64 getWidth() {
        return this->wdth;
    }

};

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
void projectOnPlane(char * ofile,
                    Plane & plane,
                    Tsph  & sph) {
    const  PS::S64 nmax = 256;
    static PS::S64 ptcl[nmax][nmax];
    static PS::F64 dens[nmax][nmax];
    static PS::F64 temp[nmax][nmax];
    static PS::F64 tmax[nmax][nmax];
    static PS::F64 xhe4[nmax][nmax];
    static PS::F64 xc12[nmax][nmax];
    static PS::F64 xo16[nmax][nmax];
    static PS::F64 xime[nmax][nmax];
    static PS::F64 xige[nmax][nmax];
    
    PS::F64 wdth  = plane.getWidth();
    PS::F64 dx    = wdth / (PS::F64)nmax;
    PS::F64 dxinv = 1. / dx;
    
    for(PS::S64 i = 0; i < nmax; i++) {
        for(PS::S64 j = 0; j < nmax; j++) {
            ptcl[i][j] = 0;
            dens[i][j] = 0.;
            temp[i][j] = 0.;
            tmax[i][j] = 0.;
            xhe4[i][j] = 0.;
            xc12[i][j] = 0.;
            xo16[i][j] = 0.;
            xime[i][j] = 0.;
            xige[i][j] = 0.;
        }
    }

    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        if(plane.calcDistance(sph[i].pos) > sph[i].ksr) {
            continue;
        }
        PS::F64vec pvec  = plane.calcCoordinateOnPlane(sph[i].pos);
        PS::S64    ix    = (PS::S64)((pvec[0] + 0.5 * wdth) * dxinv);
        PS::S64    iy    = (PS::S64)((pvec[1] + 0.5 * wdth) * dxinv);
        if(ix < 0 || nmax <= ix) {
            continue;
        }
        if(iy < 0 || nmax <= iy) {
            continue;
        }
        ptcl[ix][iy] += 1;
        dens[ix][iy] += sph[i].dens;
        temp[ix][iy] += sph[i].temp;
        tmax[ix][iy] += sph[i].tempmax[0];
        xhe4[ix][iy] += sph[i].cmps[0];
        xc12[ix][iy] += sph[i].cmps[1];
        xo16[ix][iy] += sph[i].cmps[2];
        xime[ix][iy] += (sph[i].cmps[5] + sph[i].cmps[6] + sph[i].cmps[7]
                         + sph[i].cmps[8] + sph[i].cmps[9]);
        xige[ix][iy] += (sph[i].cmps[10] + sph[i].cmps[11] + sph[i].cmps[12]);
    }

    static PS::S64 ptcl_g[nmax][nmax];
    static PS::F64 dens_g[nmax][nmax];
    static PS::F64 temp_g[nmax][nmax];
    static PS::F64 tmax_g[nmax][nmax];
    static PS::F64 xhe4_g[nmax][nmax];
    static PS::F64 xc12_g[nmax][nmax];
    static PS::F64 xo16_g[nmax][nmax];
    static PS::F64 xime_g[nmax][nmax];
    static PS::F64 xige_g[nmax][nmax];

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(ptcl, ptcl_g, nmax*nmax, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(dens, dens_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(temp, temp_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(tmax, tmax_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(xhe4, xhe4_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(xc12, xc12_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(xo16, xo16_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(xime, xime_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(xige, xige_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(PS::S64 i = 0; i < nmax; i++) {
        for(PS::S64 j = 0; j < nmax; j++) {
            PS::F64 pinv = ((ptcl_g[i][j] != 0) ? (1. / ptcl_g[i][j]) : 0.);
            dens_g[i][j] *= pinv;
            temp_g[i][j] *= pinv;
            tmax_g[i][j] *= pinv;
            xhe4_g[i][j] *= pinv;
            xc12_g[i][j] *= pinv;
            xo16_g[i][j] *= pinv;
            xime_g[i][j] *= pinv;
            xige_g[i][j] *= pinv;
        }
    }

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        for(PS::S64 i = 0; i < nmax; i++) {
            for(PS::S64 j = 0; j < nmax; j++) {
                /*
                  if(ptcl[i][j] == 0) {
                  continue;
                  }
                */
                PS::F64 px = dx * (PS::F64)i - 0.5 * wdth;
                PS::F64 pz = dx * (PS::F64)j - 0.5 * wdth;
                fprintf(fp, "%+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %8lld\n",
                        px, pz,
                        dens_g[i][j], temp_g[i][j], tmax_g[i][j],
                        xhe4_g[i][j], xc12_g[i][j], xo16_g[i][j],
                        xime_g[i][j], xige_g[i][j], ptcl_g[i][j]);
            }
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
    Plane plane;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld%lld", &ibgn, &iend);
    plane.readData(fp);
    fclose(fp);
    plane.writeData(stdout);

    for(PS::S64 itime = ibgn; itime <= iend; itime++) {        
        { // input data
            char sfile[1024];
            sprintf(sfile, "%s/sph_t%04d", idir, itime);
            sph.readParticleAscii(sfile, "%s_p%06d_i%06d.dat");
        }

        {
            char ofile[1024];
            sprintf(ofile, "%s/anim_%04d.dat", odir, itime);
            projectOnPlane(ofile, plane, sph);
        }
    }

    PS::Finalize();

    return 0;
}
