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

class Plane {
    PS::F64vec nvec;    // Normal vector of the projection surface
    PS::F64    cnst;    // Distance of the projection surface from the origin
    PS::F64vec base;    // The nearest point on the projection surface from the origin
    PS::F64vec axis[2]; // Newly defined vectors of the projection surface
    PS::F64vec pdsp;    // Vector of pallaral displacement on the projection surface
    PS::F64    wdth;    // Width of the projection surface

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
        } else if (this->nvec[0] == 1. && this->nvec[1] == 0. && this->nvec[2] == 0.) {
            this->axis[0] = PS::F64vec(0., 1., 0.);
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
        this->pdsp    = 0.;
        this->wdth    = 0.;
    }

    PS::F64 readData(FILE *fp) {
        fscanf(fp, "%lf%lf%lf%lf",
               &this->nvec[0], &this->nvec[1], &this->nvec[2], &this->cnst);
        fscanf(fp, "%lf%lf", &this->pdsp[0], &this->pdsp[1]);
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
            fprintf(fp, "pdsp: %+e %+e %+e\n",
                    this->pdsp[0], this->pdsp[1], this->pdsp[2]);
            fprintf(fp, "wdth: %+e\n", this->wdth);
        }
    }

    inline PS::F64 calcDistance(PS::F64vec pos) {
        return fabs(this->nvec * (pos - this->base));
    }

    inline PS::F64vec calcCoordinateOnPlane(PS::F64vec pos) {
        PS::F64 px = this->axis[0] * (pos - this->base) + pdsp[0];
        PS::F64 py = this->axis[1] * (pos - this->base) + pdsp[1];
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
void calcDensityCenter(Tsph  & sph,
                       PS::F64vec & xdens_g,
                       PS::F64vec & vdens_g,
                       PS::F64vec & adens_g) {
    PS::F64    tdens = 0.;
    PS::F64vec xdens = 0.;
    PS::F64vec vdens = 0.;
    PS::F64vec adens = 0.;

    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
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
void projectOnPlane(char * ofile,
                    Plane & plane,
                    Tsph  & sph) {
    const  PS::S64 nmax = 256;
    static PS::S64 ptcl[nmax][nmax];
    static PS::F64 dens[nmax][nmax];
    static PS::F64 temp[nmax][nmax];
    static PS::F64 omgz[nmax][nmax];
    //static PS::F64 domz[nmax][nmax];
    static PS::F64 istr[nmax][nmax];
    
    PS::F64 wdth  = plane.getWidth();
    PS::F64 dx    = wdth / (PS::F64)nmax;
    PS::F64 dxinv = 1. / dx;
    
    for(PS::S64 i = 0; i < nmax; i++) {
        for(PS::S64 j = 0; j < nmax; j++) {
            ptcl[i][j] = 0;
            dens[i][j] = 0.;
            temp[i][j] = 0.;
            omgz[i][j] = 0.;
            //domz[i][j] = 0.;
            istr[i][j] = 0.;
        }
    }

    PS::F64vec xdens = 0.;
    PS::F64vec vdens = 0.;
    PS::F64vec adens = 0.;
    calcDensityCenter(sph, xdens, vdens, adens);

    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].pos -= xdens;
        sph[i].vel -= vdens;
        sph[i].acc -= adens;
        PS::F64vec omg_i  = sph[i].pos ^ sph[i].vel;
        PS::F64    omgz_i = (1. / (sph[i].pos[0] * sph[i].pos[0] + sph[i].pos[1] * sph[i].pos[1])) * omg_i[2];
        //PS::F64vec trq_i  = sph[i].pos ^ sph[i].acc;
        //PS::F64    domz_i = trq_i[2] / (sph[i].pos[0] * sph[i].pos[0] + sph[i].pos[1] * sph[i].pos[1]);
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
        omgz[ix][iy] += omgz_i;
        //domz[ix][iy] += domz_i;
        istr[ix][iy] += (sph[i].istar + 1);
    }

    static PS::S64 ptcl_g[nmax][nmax];
    static PS::F64 dens_g[nmax][nmax];
    static PS::F64 temp_g[nmax][nmax];
    static PS::F64 omgz_g[nmax][nmax];
    //static PS::F64 domz_g[nmax][nmax];
    static PS::F64 istr_g[nmax][nmax];

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(ptcl, ptcl_g, nmax*nmax, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(dens, dens_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(temp, temp_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(omgz, omgz_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //ierr = MPI_Allreduce(domz, domz_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(istr, istr_g, nmax*nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(PS::S64 i = 0; i < nmax; i++) {
        for(PS::S64 j = 0; j < nmax; j++) {
            PS::F64 pinv = ((ptcl_g[i][j] != 0) ? (1. / ptcl_g[i][j]) : 0.);
            dens_g[i][j] *= pinv;
            temp_g[i][j] *= pinv;
            omgz_g[i][j] *= pinv;
            //domz_g[i][j] *= pinv;
            istr_g[i][j] *= pinv;
        }
    }

    if(PS::Comm::getRank() == 0) {
        FILE * fp = fopen(ofile, "w");
        for(PS::S64 i = 0; i < nmax; i++) {
            for(PS::S64 j = 0; j < nmax; j++) {
                PS::F64 px = dx * (PS::F64)i - 0.5 * wdth;
                PS::F64 pz = dx * (PS::F64)j - 0.5 * wdth;
#if 0
                fprintf(fp, "%+e %+e %+e %+e %+e %+e %+e\n",
                        px, pz,
                        dens_g[i][j], temp_g[i][j], omgz_g[i][j],
                        domz_g[i][j], istr_g[i][j]);
#endif
                fprintf(fp, "%+e %+e %+e %+e %+e %+e\n",
                        px, pz,
                        dens_g[i][j], temp_g[i][j], omgz_g[i][j],
                        istr_g[i][j]);
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
    PS::S64 ibgn, iend, dsnp;
    Plane plane;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lld%lld%lld", &ibgn, &iend, &dsnp);
    plane.readData(fp);
    fclose(fp);
    plane.writeData(stdout);

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
        sprintf(ofile, "%s/anim_%04d.dat", odir, itime);
        projectOnPlane(ofile, plane, sph);
    }

    PS::Finalize();

    return 0;
}
