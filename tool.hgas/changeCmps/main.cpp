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
#include "vector_x86.hpp"
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_hgas.hpp"

class SPH3D : public HelmholtzGas {
public:
    SPH3D() {
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

    void read(FILE * fp) {
        fscanf(fp, "%lld%lld%lf", &this->id, &this->istar, &this->mass);         //  3
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
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) { // 32 -- 44
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        /*
        fscanf(fp, "%d", &this->fnse);
        */
        fscanf(fp, "%lf%lf%lf", &this->tempmax[0], &this->tempmax[1], &this->tempmax[2]);
    }

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
};

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    PS::S64 nptcl;
    char ifile[64];
    char ofile[64];
    PS::S64 flag;
    PS::F64 rmin;
    NR::Nucleon cmps;

    {
        FILE * fp = fopen(argv[1], "r");
        fscanf(fp, "%lld%s%s", &nptcl, ifile, ofile);
        fscanf(fp, "%lf", &rmin);
        fscanf(fp, "%lf%lf%lf%lf%lf", &cmps[0], &cmps[1], &cmps[2], &cmps[3], &cmps[4]);
        printf("nptcl: %8d\n", nptcl);
        printf("ifile: %s\n", ifile);
        printf("ofile: %s\n", ofile);
        printf("rmin: %+e\n", rmin);
        printf("He: %+.3e C: %+.3e O: %+.3e Ne: %+.3e Mg: %+.3e\n",
               cmps[0], cmps[1], cmps[2], cmps[3], cmps[4]);
        fclose(fp);
    }    

    {
        FILE * ifp = fopen(ifile, "r");
        assert(ifp);
        FILE * ofp = fopen(ofile, "w");
        for(PS::S64 i = 0; i < nptcl; i++) {
            SPH3D sph;
            sph.read(ifp);
            PS::F64 r2 = sph.pos * sph.pos;
            if(r2 >= rmin * rmin) {
                sph.cmps = cmps;
            }
            sph.write(ofp);
        }
        fclose(ifp);
        fclose(ofp);
    }

    MPI_Finalize();

    return 0;
}
