#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>
#include <mpi.h>
#include <Eigen/Eigenvalues>

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

    /*
    void write(FILE * fp = stdout) {
        fprintf(fp, "%8d %2d %+e",  this->id, this->istar, this->mass);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->vel[0], this->vel[1], this->vel[2]);
        fprintf(fp, " %+e %+e %+e", this->uene, this->alph, this->alphu);
        fprintf(fp, " %+e", this->ksr);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        fprintf(fp, "\n");
    }
    */

    inline PS::F64mat calcMomentOfInertia(PS::F64vec xc) {
        PS::F64mat qq = 0.;
        PS::F64vec dx = this->pos - xc;
        PS::F64 r2 = dx * dx;
        qq.xx = this->mass * (r2 - dx[0] * dx[0]);
        qq.xy = this->mass * (   - dx[0] * dx[1]);
        qq.xz = this->mass * (   - dx[0] * dx[2]);
        qq.yy = this->mass * (r2 - dx[1] * dx[1]);
        qq.yz = this->mass * (   - dx[1] * dx[2]);
        qq.zz = this->mass * (r2 - dx[2] * dx[2]);
        return qq;
    }
};


class BlackHoleNeutronStarAnalysis : public BlackHoleNeutronStar {
public:
    BlackHoleNeutronStarAnalysis() {
        this->id    = 0;
        this->istar = 0;
        this->mass  = 0.;
        this->pos   = 0.;
        this->vel   = 0.;
        this->acc   = 0.;
        this->eps   = 0.;
        this->pot   = 0.;
    }

    void readAscii(FILE * fp) {
        fscanf(fp, "%lld%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->istar, &this->mass,
               &this->pos[0], &this->pos[1], &this->pos[2],
               &this->vel[0], &this->vel[1], &this->vel[2],
               &this->acc[0], &this->acc[1], &this->acc[2],
               &this->eps, &this->pot);
    }

    inline PS::F64mat calcMomentOfInertia(PS::F64vec xc) {
        PS::F64mat qq = 0.;
        PS::F64vec dx = this->pos - xc;
        PS::F64 r2 = dx * dx;
        qq.xx = this->mass * (r2 - dx[0] * dx[0]);
        qq.xy = this->mass * (   - dx[0] * dx[1]);
        qq.xz = this->mass * (   - dx[0] * dx[2]);
        qq.yy = this->mass * (r2 - dx[1] * dx[1]);
        qq.yz = this->mass * (   - dx[1] * dx[2]);
        qq.zz = this->mass * (r2 - dx[2] * dx[2]);
        return qq;
    }
};

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);
    PS::ParticleSystem<BlackHoleNeutronStarAnalysis> bhns;
    bhns.initialize();
    bhns.createParticle(0);
    bhns.setNumberOfParticleLocal(0);

    char idir[1024], odir[1024];
    PS::F64 rcrit;
    PS::S64 ibgn, iend;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lf", &rcrit);
    fscanf(fp, "%lld%lld", &ibgn, &iend);
    fclose(fp);

    char ofile[1024];
    sprintf(ofile, "%s/axis.dat", odir);
    FILE * fout = fopen(ofile, "w");
    fprintf(fout, "# idir: %s\n", idir);
    fprintf(fout, "# crit: %+e\n", rcrit);
    fprintf(fout, "# time, Ixx, Iyy, Izz, Ixx, Iyy, Izz\n");
    for(PS::S64 itime = ibgn; itime <= iend; itime++) {
        char sfile[1024], bfile[1024];
        sprintf(sfile, "%s/sph_t%04d.dat",  idir, itime);
        sprintf(bfile, "%s/bhns_t%04d.dat", idir, itime);
        fp = fopen(sfile, "r");
        assert(fp);
        sph.readParticleAscii(sfile);
        fclose(fp);
        fp = fopen(bfile, "r");
        assert(fp);
        bhns.readParticleAscii(bfile);
        fclose(fp);

        PS::F64mat qq[2] = {0., 0.};
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            PS::F64mat qtemp = sph[i].calcMomentOfInertia(bhns[0].pos);
            qq[0] = qq[0] + qtemp;
            PS::F64vec dx = sph[i].pos - bhns[0].pos;
            PS::F64    r2 = dx * dx;
            if(r2 < rcrit * rcrit) {
                qq[1] = qq[1] + qtemp;
            }
        }

        fprintf(fout, "%5d", itime);
        Eigen::Matrix3d MomentOfInertia[2];
        for(PS::S64 k = 0; k < 2; k++) {
            MomentOfInertia[k] <<
                qq[k].xx, qq[k].xy, qq[k].xz,
                qq[k].xy, qq[k].yy, qq[k].yz,
                qq[k].xz, qq[k].yz, qq[k].zz;
            Eigen::EigenSolver<Eigen::Matrix3d> es(MomentOfInertia[k], false);
            std::vector<PS::F64> ev;
            ev.push_back(es.eigenvalues()[0].real());
            ev.push_back(es.eigenvalues()[1].real());
            ev.push_back(es.eigenvalues()[2].real());
            std::sort(ev.begin(), ev.end());
            fprintf(fout, " %+e %+e %+e", ev[0], ev[1], ev[2]);
        }
        fprintf(fout, "\n");

    }
    fclose(fout);

    MPI_Finalize();

    return 0;
}
