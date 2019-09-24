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

MPI_Datatype ParticleRecordMPI;

class ParticleRecord {
public:
    ParticleRecord() {
        this->id   = 0;
        this->time = 0.;
        this->temp = 0.;
        this->dens = 0.;
    }

    PS::S64 id;
    PS::F64 time;
    PS::F64 temp;
    PS::F64 dens;

    bool readLine(FILE * fp) {
        PS::S64 ret = fscanf(fp, "%lf%lld%lf%lf", &this->time, &this->id, &this->temp, &this->dens);
        return ((ret == 4) ? true : false);
    }

    void writeLine(FILE * fp) {
        fprintf(fp, "%16.10f %10d %+e %+e\n", this->time, this->id, this->temp, this->dens);
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

void createMPIDataTypeOfParticleRecord() {
    PS::S32 ret;
    PS::S32 block_lengths[4] = {1, 1, 1, 1};
    MPI_Aint displacement[4] = {0, 8, 16, 24};
    MPI_Datatype typelist[4] = {MPI_LONG_LONG_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    ret = MPI_Type_struct(4, block_lengths, displacement, typelist, &ParticleRecordMPI);
    assert(ret == MPI_SUCCESS);
    ret = MPI_Type_commit(&ParticleRecordMPI);
    assert(ret == MPI_SUCCESS);
}

template <class TypeSort>
void sortQuicklyParticleRecord(TypeSort * prec,
                               PS::S64 left,
                               PS::S64 right) {
    PS::S64 i    = left;
    PS::S64 j    = right;
    PS::S64 axis = prec[(i+j)/2].id;
    while(1) {
        while(prec[i].id < axis) i++;
        while(prec[j].id > axis) j--;
        if(i >= j) break;
        TypeSort temp = prec[i];
        prec[i] = prec[j];
        prec[j] = temp;
        i++;
        j--;
    }
    if(left < i-1) {
        sortQuicklyParticleRecord(prec, left, i-1);
    }
    if(j+1 < right) {
        sortQuicklyParticleRecord(prec, j+1, right);
    }
}

template <class TypeSort>
void sortQuicklyParticleRecordTime(TypeSort * prec,
                                   PS::S64 left,
                                   PS::S64 right) {
    PS::S64 i    = left;
    PS::S64 j    = right;
    PS::F64 axis = prec[(i+j)/2].time;
    while(1) {
        while(prec[i].time < axis) i++;
        while(prec[j].time > axis) j--;
        if(i >= j) break;
        TypeSort temp = prec[i];
        prec[i] = prec[j];
        prec[j] = temp;
        i++;
        j--;
    }
    if(left < i-1) {
        sortQuicklyParticleRecordTime(prec, left, i-1);
    }
    if(j+1 < right) {
        sortQuicklyParticleRecordTime(prec, j+1, right);
    }
}

PS::S64 getRankToBeDistributed(PS::S64 id,
                               PS::S64 id1st,
                               PS::S64 idint,
                               PS::S64 nptcl) {
    PS::S64 add   = (id - id1st) / idint;
    PS::S64 iproc = add * PS::Comm::getNumberOfProc() / nptcl;
    return iproc;
}

int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);
    createMPIDataTypeOfParticleRecord();

    const PS::S32 nprocmax = 16384;
    static PS::S32 sendcount[nprocmax];
    static PS::S32 senddspls[nprocmax];
    static PS::S32 recvcount[nprocmax];
    static PS::S32 recvdspls[nprocmax];
    assert(PS::Comm::getNumberOfProc() <= nprocmax);

    char idir[1024], ifile[1024];
    char odir[1024], ofile[1024];
    PS::S64 nptcl;
    PS::S64 id1st, idint;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s%s%lld", idir, odir, &nptcl);
    fscanf(fp, "%lld%lld", &id1st, &idint);
    fclose(fp);

    // Count line in advance
    PS::S64 nline = 0;
    for(PS::S64 itime = 0; itime <= 99; itime++) {
        char tfile[1024];
        FILE *fp = NULL;
        PS::S64 tdir = 0;
        sprintf(tfile, "%s/t%02d/debug_p%06d.log", idir, itime, PS::Comm::getRank());
        fp = fopen(tfile, "r");
        if(fp == NULL) {
            if(PS::Comm::getRank() == 0) fprintf(stderr, "Not found %s\n", tfile);
            continue;
        }
        ParticleRecord prectmp;
        while(prectmp.readLine(fp)) nline++;
        fclose(fp);
    }
    PS::S64 nlinetot = PS::Comm::getSum(nline);
    PS::S64 nstep    = nlinetot / nptcl;

    // Allocate memory
    PS::S64 nlineloc = ((nline / 64) + 1) * 64;
    PS::S64 nlineglb = ((nptcl / PS::Comm::getNumberOfProc() + 1) * nstep / 64 + 1) * 64;
    ParticleRecord * precloc;
    ParticleRecord * precglb;
    int ret;
    ret = posix_memalign((void **)&precloc, 64, sizeof(ParticleRecord)*nlineloc);
    assert(ret == 0);
    ret = posix_memalign((void **)&precglb, 64, sizeof(ParticleRecord)*nlineglb);
    assert(ret == 0);

    // Read data
    for(PS::S64 itime = 0, i = 0; itime <= 99; itime++) {
        char tfile[1024];
        FILE *fp = NULL;
        PS::S64 tdir = 0;
        sprintf(tfile, "%s/t%02d/debug_p%06d.log", idir, itime, PS::Comm::getRank());
        fp = fopen(tfile, "r");
        if(fp == NULL) {
            continue;
        }
        while(precloc[i].readLine(fp)) i++;
        fclose(fp);
    }

    // Sort ID
    if(nline > 0) {
        sortQuicklyParticleRecord(precloc, 0, nline-1);
    }

    // Prepare data distribution
    for(PS::S64 i = 0; i < nline; i++) {
        PS::S64 iproc = getRankToBeDistributed(precloc[i].id, id1st, idint, nptcl);
        sendcount[iproc]++;
    }
    for(PS::S64 i = 0, nsend = 0; i < PS::Comm::getNumberOfProc()+1; i++) {
        senddspls[i] = nsend;
        nsend += sendcount[i];
    }
    ret = MPI_Alltoall(sendcount, 1, MPI_INT,
                       recvcount, 1, MPI_INT,
                       MPI_COMM_WORLD);
    assert(ret == MPI_SUCCESS);
    for(PS::S64 i = 0, nrecv = 0; i < PS::Comm::getNumberOfProc()+1; i++) {
        recvdspls[i] = nrecv;
        nrecv += recvcount[i];
    }

    // Distribute data
    ret = MPI_Alltoallv(precloc, sendcount, senddspls, ParticleRecordMPI,
                        precglb, recvcount, recvdspls, ParticleRecordMPI,
                        MPI_COMM_WORLD);
    assert(ret == MPI_SUCCESS);

    // Sort time for each particle
    PS::S64 nrecv = recvdspls[PS::Comm::getNumberOfProc()];
    if(nrecv > 0) {
        sortQuicklyParticleRecord(precglb, 0, nrecv-1);
    }
    for(PS::S64 i = 0; i < nrecv; i += nstep) {
        sortQuicklyParticleRecordTime(precglb, i, i+nstep-1);
    }

    // Output particle data for each file
    sprintf(ofile, "%s/prec%06d.dat", odir, PS::Comm::getRank());
    fp = fopen(ofile, "w");
    PS::S64 pid   = -1;
    PS::F64 ptime = -1.;
    for(PS::S64 i = 0; i < nrecv; i++) {
        if(precglb[i].id == pid && precglb[i].time == ptime) {
            continue;
        }
        precglb[i].writeLine(fp);
        pid   = precglb[i].id;
        ptime = precglb[i].time;
    }
    fclose(fp);

    free(precloc);
    free(precglb);
    
    PS::Finalize();

    return 0;
}
