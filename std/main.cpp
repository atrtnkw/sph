#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};

#include "particle_simulator.hpp"
#include "hdr_time.hpp"
#include "hdr_run.hpp"
#include "vector_x86.hpp"
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_igas.hpp"
#include "hdr_density.hpp"
#include "hdr_hydro.hpp"

template <class Tsph>
void referEquationOfState(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].referEquationOfState();
    }
}

template <class Tsph>
void calcBalsaraSwitch(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].calcBalsaraSwitch();
    }
}

template <class Tsph>
void calcAlphaDot(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].calcAlphaDot();
    }
}

template <class Tsph>
PS::F64 calcEnergy(Tsph & sph) {
    PS::F64 eloc = 0.;
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        eloc += sph[i].calcEnergy();
    }
    PS::F64 eglb = PS::Comm::getSum(eloc);
    return eglb;
}

template <class Tsph>
PS::F64vec calcMomentum(Tsph & sph) {
    PS::F64vec momloc = 0.;
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        momloc += sph[i].calcMomentum();
    }
    PS::F64vec momglb = PS::Comm::getSum(momloc);
    return momglb;
}

template <class Tsph>
PS::F64 calcTimeStep(Tsph & sph) {
    PS::S32 nloc = sph.getNumberOfParticleLocal();
    PS::F64 dtc = RP::MaximumTimestep;
    for(PS::S32 i = 0; i < nloc; i++) {
        PS::F64 dttmp = sph[i].calcTimeStep();
        if(dttmp < dtc)
            dtc = dttmp;
    }
    dtc = PS::Comm::getMinValue(dtc);

    PS::F64 dt = RP::Timestep;
    if(dt > dtc) {
        while(dt > dtc) {
            dt *= 0.5;
            if(dt <= RP::MinimumTimestep) {
                printf("Time: %.10f Timestep: %+e\n", RP::Time, dtc);
                sph.writeParticleAscii("snap/smalldt.dat");
                PS::Finalize();
                exit(0);
            }
        }        
    } else if(2. * dt <= dtc && 2. * dt <= RP::MaximumTimestep) {
        PS::F64 dt2 = 2. * dt;
        if(RP::Time - (PS::S64)(RP::Time / dt2) * dt2 == 0.) {
            dt *= 2.;
        }
    }

    return dt;
}

template <class Tdinfo,
          class Tsph,
          class Tdensity,
          class Thydro>
void calcSPHKernel(Tdinfo & dinfo,
                   Tsph & sph,
                   Tdensity & density,
                   Thydro & hydro) {
    calcDensityKernel(dinfo, sph, density);
    referEquationOfState(sph);
    calcBalsaraSwitch(sph);
    hydro.calcForceAllAndWriteBack(calcHydro(), sph, dinfo);
    calcAlphaDot(sph);
}

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::DomainInfo dinfo;
    dinfo.initialize();
    PS::ParticleSystem<GeneralSPH> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);
    PS::TreeForForceShort<Density, DensityEPI, DensityEPJ>::Gather density;
    density.initialize(0);
    PS::TreeForForceShort<Hydro, HydroEPI, HydroEPJ>::Symmetry hydro;
    hydro.initialize(0);

    initializeSimulation();

    if(atoi(argv[1]) == 0) {
        startSimulation(argv, dinfo, sph, density, hydro);
    } else {
        restartSimulation(argv, dinfo, sph);
    }

    loopSimulation(dinfo, sph, density, hydro);

    finalizeSimulation(sph);

    PS::Finalize();

    return 0;
}
