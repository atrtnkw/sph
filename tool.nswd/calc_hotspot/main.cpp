#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>
//#include <mpi.h>
//#include <Eigen/Eigenvalues>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};
class Physics;

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

    void copyFromForce(const Physics & physics);
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

class MeshData : public SPHAnalysis {
public:
    PS::F64 tempu;

    MeshData() {
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

    void writeAscii(FILE *fp) {
        fprintf(fp, "%8d", this->id);
        fprintf(fp, " %+e %+e %+e", this->pos[0], this->pos[1], this->pos[2]);
        fprintf(fp, " %+e %+e %+e", this->dens, this->temp, this->tempu);
        fprintf(fp, "\n");
    }

    void copyFromForce(const Physics & physics);
};

template <class Tmesh>
void generateMeshData(Tmesh & mesh) {
    if(PS::Comm::getRank() != 0) {
        return;
    }
    PS::F64 xmin = -5e8;
    PS::F64 xmax = +5e8;
    PS::S64 nx   = 100;
    PS::S64 nx3  = nx * nx * nx;
    PS::F64 dx   = (xmax - xmin) / (PS::F64)nx;
    mesh.setNumberOfParticleLocal(nx3);
    for(PS::S64 i = 0; i < nx; i++) {
        for(PS::S64 j = 0; j < nx; j++) {
            for(PS::S64 k = 0; k < nx; k++) {
                PS::S64 id = i * nx * nx + j * nx + k;
                mesh[id].id     = id;
                mesh[id].mass   = 0.;
                mesh[id].pos[0] = xmin + dx * (PS::F64)i;
                mesh[id].pos[1] = xmin + dx * (PS::F64)j;
                mesh[id].pos[2] = xmin + dx * (PS::F64)k;
                mesh[id].ksr    = 1e7;
                mesh[id].dens   = 0.;
                mesh[id].temp   = 0.;
            }
        }
    }
}

class Physics {
public:
    PS::F64 dens;
    PS::F64 temp;
    PS::F64 uene;
    void clear() {
        dens = 0.;
        temp = 0.;
        uene = 0.;
    }
};

void SPHAnalysis::copyFromForce(const Physics & physics) {
    this->dens = physics.dens;
    this->temp = physics.temp;
    this->uene = physics.uene;
}

void MeshData::copyFromForce(const Physics & physics) {
    this->dens = physics.dens;
    this->temp = physics.temp;
    this->uene = physics.uene;
}

class PhysicsEPI {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    bool       meshornot;

    void copyFromFP(const SPHAnalysis & sph) {
        this->id   = sph.id;
        this->mass = sph.mass;        
        this->pos  = sph.pos;
        this->meshornot = false;
    }

    void copyFromFP(const MeshData & mesh) {
        this->id   = mesh.id;
        this->mass = mesh.mass;        
        this->pos  = mesh.pos;
        this->meshornot = true;
    }

    PS::F64vec getPos() const {
        return this->pos;
    }

};

class PhysicsEPJ {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    dinv;
    PS::F64    temp;
    PS::F64    uene;
    bool       meshornot;

    void copyFromFP(const SPHAnalysis & sph) {
        this->id   = sph.id;
        this->mass = sph.mass;        
        this->pos  = sph.pos;
        this->ksr  = sph.ksr;
        this->hinv = 1. / sph.ksr;
        this->dinv = 1. / sph.dens;
        this->temp = sph.temp;
        this->uene = sph.uene;
        this->meshornot = false;
    }

    void copyFromFP(const MeshData & mesh) {
        this->id   = mesh.id;
        this->mass = mesh.mass;        
        this->pos  = mesh.pos;
        this->ksr  = mesh.ksr;
        this->hinv = 1. / mesh.ksr;
        this->dinv = 1. / mesh.dens;
        this->temp = mesh.temp;
        this->uene = mesh.uene;
        this->meshornot = true;
    }

    PS::F64vec getPos() const {
        return this->pos;
    }

    PS::F64 getRSearch() const {
        return this->ksr;
    }

    void setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
};

struct calcPhysics {

    const PS::F64 ceff0 = +3.342253804929802286e+00;
    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = ((1. - r > 0.) ? (1. - r) : 0.);
        PS::F64 rmin2 = rmin * rmin;
        return ceff0 * rmin2 * rmin2 * (1. + 4. * r);
    }

    void operator () (const PhysicsEPI * epi,
                      const PS::S32 nip,
                      const PhysicsEPJ * epj,
                      const PS::S32 njp,
                      Physics * physics) {

        for(PS::S64 i = 0; i < nip; i++) {

            if(!epi[i].meshornot) {
                continue;
            }

            PS::F64 dens = 0.;
            PS::F64 temp = 0.;
            PS::F64 uene = 0.;

            for(PS::S64 j = 0; j < njp; j++) {
                PS::F64vec dx  = epi[i].pos - epj[j].pos;
                PS::F64    r2  = dx * dx;
                PS::F64    r1  = sqrt(r2);
                PS::F64    q   = r1 * epj[j].hinv;
                PS::F64    kw0 = kernel0th(q);
                PS::F64    hi3 = epj[j].hinv * epj[j].hinv * epj[j].hinv;

                PS::F64 densi = epj[j].mass * hi3 * kw0;
                PS::F64 tempi = epj[j].mass * hi3 * kw0 * epj[j].temp * epj[j].dinv;
                PS::F64 uenei = epj[j].mass * hi3 * kw0 * epj[j].uene * epj[j].dinv;
                dens += ((epj[j].mass != 0.) ? densi : 0.);
                temp += ((epj[j].mass != 0.) ? tempi : 0.);
                uene += ((epj[j].mass != 0.) ? uenei : 0.);
            }

            physics[i].dens = dens;
            physics[i].temp = temp;
            physics[i].uene = uene;
        }
        
    }

};

template <class Tdinfo,
          class Tsph,
          class Tmesh,
          class Tphysics>
void calcMeshData(Tdinfo & dinfo,
                  Tsph & sph,
                  Tmesh & mesh,
                  Tphysics & physics) {
    PS::MT::init_genrand(0);
    dinfo.decomposeDomainAll(sph);
    sph.exchangeParticle(dinfo);
    mesh.exchangeParticle(dinfo);
    physics.setParticleLocalTree(mesh);
    physics.setParticleLocalTree(sph, false);
    physics.calcForceMakingTree(calcPhysics(), dinfo);
    PS::S64 nmesh = mesh.getNumberOfParticleLocal();
    for(PS::S64 i = 0; i < nmesh; i++) {
        mesh[i].copyFromForce(physics.getForce(i));
    }
    for(PS::S64 i = 0; i < nmesh; i++) {
        NR::Nucleon cmps;
        cmps[1] = cmps[2] = 0.5;
        PS::F64 dens = mesh[i].dens;
        PS::F64 uene = mesh[i].uene;
        PS::F64 temp = mesh[i].temp;
        PS::F64 tout = 0.;
        PS::F64 pout, cout, sout;
        if(dens != 0.) {
            flash_helmholtz_(&dens, &uene, &temp, cmps.getPointer(),
                             &pout, &cout, &tout, &sout);
        }
        mesh[i].tempu = tout;
    }
}

template <class Tsph,
          class Tbhns>
void shiftCenter(Tsph & sph,
                 Tbhns & bhns) {
    PS::S64    nbhns = bhns.getNumberOfParticleLocal();
    PS::F64vec ploc  = 0.;
    if(nbhns == 1) {
        ploc = bhns[0].pos;
    } else {
        ploc = 0.;
    }
    PS::F64vec pbhns = PS::Comm::getSum(ploc);

    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].pos -= pbhns;
    }
}

int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);

    init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);

    PS::DomainInfo dinfo;
    dinfo.initialize();
    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);
    PS::ParticleSystem<BlackHoleNeutronStarAnalysis> bhns;
    bhns.initialize();
    bhns.createParticle(0);
    bhns.setNumberOfParticleLocal(0);
    PS::ParticleSystem<MeshData> mesh;
    mesh.initialize();
    mesh.createParticle(0);
    mesh.setNumberOfParticleLocal(0);
    PS::TreeForForceShort<Physics, PhysicsEPI, PhysicsEPJ>::Scatter physics;
    physics.initialize(0);

    if(PS::Comm::getRank() == 0) {
        printf("AT comments: NPROCESS %8d\n", PS::Comm::getNumberOfProc());
        printf("AT comments: NTHREAD  %8d\n", PS::Comm::getNumberOfThread());
    }

    char idir[1024], otype[1024];
    PS::S64 ibgn, iend;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", otype);
    fscanf(fp, "%lld%lld", &ibgn, &iend);
    fclose(fp);

    generateMeshData(mesh);

    for(PS::S64 itime = ibgn; itime <= iend; itime++) {
        char sfile[1024], bfile[1024];
        sprintf(sfile, "%s/sph_t%04d.dat",  idir, itime);
        sprintf(bfile, "%s/bhns_t%04d.dat", idir, itime);
        fp = fopen(sfile, "r");
        if(fp == NULL) {
            continue;
        }
        sph.readParticleAscii(sfile);
        fclose(fp);
        fp = fopen(bfile, "r");
        if(fp == NULL) {
            continue;
        }
        bhns.readParticleAscii(bfile);
        fclose(fp);

        shiftCenter(sph, bhns);
        calcMeshData(dinfo, sph, mesh, physics);

        char ofile[1024];
        sprintf(ofile, "%s_t%04d.dat", otype, itime);
        mesh.writeParticleAscii(ofile);
    }

    PS::Finalize();

    return 0;
}
