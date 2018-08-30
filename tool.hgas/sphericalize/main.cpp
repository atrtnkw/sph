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

    PS::F64 calcRadialVelocity() {
        PS::F64 r2 = this->pos * this->pos;
        PS::F64 rv = this->pos * this->vel;
        return (rv / sqrt(r2));
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

class Shell {
public:
    PS::S64 nump; // The number of particles in the shell
    PS::F64 mass; // The total mass in the shell
    PS::F64 ener; // The total internal energy in the shell
    PS::F64 rvel; // Radial velocity in the shell
    NR::Nucleon mele;

    Shell() {
        this->nump = 0;
        this->mass = 0.;
        this->ener = 0.;
        this->rvel = 0.;
        for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
            this->mele[k] = 0.;
        }
    }

    void increment(SPHAnalysis & sph) {
        this->nump += 1;
        this->mass += sph.mass;
        this->ener += sph.mass * sph.uene;
        this->rvel += sph.calcRadialVelocity();
        for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
            this->mele[k] += sph.mass * sph.cmps[k];
        }
    }

    void reduceAddition(Shell sloc) {
        this->nump = PS::Comm::getSum(sloc.nump);
        this->mass = PS::Comm::getSum(sloc.mass);
        this->ener = PS::Comm::getSum(sloc.ener);
        this->rvel = PS::Comm::getSum(sloc.rvel);
        for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
            this->mele[k] = PS::Comm::getSum(sloc.mele[k]);
        }
    }

};

template <class Tsph>
void sphericalizeWhiteDwarf(PS::F64 rmax,
                            PS::F64 drad,
                            char * ofile,
                            Tsph & sph) {

    PS::F64vec cntr_loc = 0.;
    PS::F64vec cvel_loc = 0.;
    PS::F64    dens_loc = 0.;
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        dens_loc += sph[i].dens;
        cntr_loc += sph[i].dens * sph[i].pos;
        cvel_loc += sph[i].dens * sph[i].vel;
    }
    PS::F64vec cntr_glb = PS::Comm::getSum(cntr_loc);
    PS::F64vec cvel_glb = PS::Comm::getSum(cvel_loc);
    PS::F64    dens_glb = PS::Comm::getSum(dens_loc);
    cntr_glb = (1. / dens_glb) * cntr_glb;
    cvel_glb = (1. / dens_glb) * cvel_glb;
    
    PS::S64 nbin = (PS::S64)(rmax / drad) + 1;
    Shell * sloc = (Shell *)malloc(sizeof(Shell) * nbin);
    
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].pos = sph[i].pos - cntr_glb;
        sph[i].vel = sph[i].vel - cvel_glb;
        PS::F64 rad1  = sqrt(sph[i].pos * sph[i].pos);
        PS::S64 ibin = (PS::S64)(rad1 / drad);
        if(ibin >= nbin) {
            continue;
        }
        ///////////// Ad hoc method //////////////
        if(sph[i].istar == 0 && sph[i].cmps[1] > 0.2 && rad1 < 3e10) {
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                if(k != 12) {
                    sph[i].cmps[k] = 0.;
                } else {
                    sph[i].cmps[k] = 1.;
                }
            }
            //continue;
        }
        //////////////////////////////////////////
        sloc[ibin].increment(sph[i]);
    }

    Shell * sglb = (Shell *)malloc(sizeof(Shell) * nbin);
    for(PS::S64 ibin = 0; ibin < nbin; ibin++) {
        sglb[ibin].reduceAddition(sloc[ibin]);
    }
    
    if(PS::Comm::getRank() == 0) {
        PS::S64 nshl = 0;
        PS::F64 mass = 0.;
        PS::F64 ener = 0.;
        PS::F64 rvel = 0.;
        PS::F64 rad0 = 0.;
        PS::F64 menc = 0.;
        NR::Nucleon mele;
        FILE * fp = fopen(ofile, "w");
        for(PS::S64 ibin = 0; ibin < nbin; ibin++) {
            nshl += sglb[ibin].nump;
            mass += sglb[ibin].mass;
            ener += sglb[ibin].ener;
            rvel += sglb[ibin].rvel;
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                mele[k] += sglb[ibin].mele[k];
            }
            if(nshl < 10000 && ibin != nbin - 1) {
            //if(nshl < 1000 && ibin != nbin - 1) {
                continue;
            }

            PS::F64 rad1   = drad * (ibin + 1);
            PS::F64 radius = ((rad0 != 0.) ? sqrt(rad0 * rad1) : 0.5 * rad1);
            PS::F64 volume = 4. * M_PI / 3. * (rad1 * rad1 * rad1 - rad0 * rad0 * rad0);
            PS::F64 dens   = mass / volume;
            PS::F64 uene   = ener / mass;
            rvel  = rvel / (PS::S64)nshl;
            NR::Nucleon cmps;
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                cmps[k] = mele[k] / mass;
            }
            menc += mass;

            PS::F64 tin = 1e9;
            PS::F64 pout, cout, tout, sout, uout;
            flash_helmholtz_(&dens, &uene, &tin, cmps.getPointer(), &pout, &cout, &tout, &sout);
            tout = (tout < 1e7) ? 1e7 : tout;
            flash_helmholtz_e_(&dens, &tout, cmps.getPointer(), &uout);
            flash_helmholtz_(&dens, &uout, &tout, cmps.getPointer(), &pout, &cout, &tout, &sout);
            fprintf(fp, "%+e %+e %+e %+e %+e %+e %+e", 
                    rad1, dens, pout, tout, uout, menc/CodeUnit::SolarMass, rvel);
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                fprintf(fp, " %+e", cmps[k]);
            }
            fprintf(fp, " %+e", sout);
            fprintf(fp, "\n");
            
            nshl = 0;
            mass = 0.;
            ener = 0.;
            rvel = 0.;
            for(PS::S64 k = 0; k < NR::NumberOfNucleon; k++) {
                mele[k] = 0.;
            }
            rad0 = rad1;
        }
        fclose(fp);
    }

    free(sloc);
    free(sglb);
                
}

int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);

    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);

    char idir[1024], odir[1024];
    PS::F64 rmax, drad;
    PS::S64 ibgn, iend;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", odir);
    fscanf(fp, "%lf%lf", &rmax, &drad);
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
            fprintf(stderr, "%s is not found.\n", tfile);
            continue;
        }
        fclose(fp);
        
        char sfile[1024];
        sprintf(sfile, "%s/t%02d/sph_t%04d", idir, tdir, itime);
        sph.readParticleAscii(sfile, "%s_p%06d_i%06d.dat");

        char ofile[1024];
        sprintf(ofile, "%s/poly.dat", odir);
        sphericalizeWhiteDwarf(rmax, drad, ofile, sph);

    }

    PS::Finalize();

    return 0;
}
