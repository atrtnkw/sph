#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};
inline void debugFunction(const char * const literal);
template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls,
          class Tdensity,
          class Thydro,
          class Tgravity>
void calcSPHKernel(Tdinfo & dinfo,
                   Tsph & sph,
                   Tbhns & bhns,
                   Tmsls & msls,
                   Tdensity & density,
                   Thydro & hydro,
                   Tgravity & gravity);

#include "particle_simulator.hpp"
#include "hdr_time.hpp"
#include "hdr_run.hpp"
#include "vector_x86.hpp"
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#ifdef USE_IDEAL
#include "hdr_igas.hpp"
#elif defined USE_HELMHOLTZ
#include "hdr_hgas.hpp"
#else
#error Use what eos is ?
#endif
#include "hdr_msls.hpp"
#include "hdr_bhns.hpp"
#include "hdr_density.hpp"
#include "hdr_hydro.hpp"
#include "hdr_gravity.hpp"

inline void debugFunction(const char * const literal) {
    if(PS::Comm::getRank() == 0) {
        printf("hogehogehoge %.16e %s\n", RP::Time, literal);
        fflush(stdout);
    }
}

template <class Tsph>
void calcAbarZbar(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].calcAbarZbar();
    }
}

template <class Tsph>
void referEquationOfState(Tsph & sph) {
    if(RP::FlagDamping != 1) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].referEquationOfState();
        }
    } else {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].referEquationOfStateDamping1();
        }
    }
}

template <class Tsph>
void calcReleasedNuclearEnergy(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].calcReleasedNuclearEnergy();        
    }
}

/*
template <class Tsph>
void predictReleasedNuclearEnergy(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].dens0 = sph[i].dens;
        sph[i].temp0 = sph[i].temp;
        sph[i].cmps0 = sph[i].cmps;
        sph[i].dnuc = GeneralSPH::calcReleasedNuclearEnergy(RP::Timestep,
                                                            sph[i].dens,
                                                            sph[i].temp,
                                                            sph[i].cmps.getPointer());
        //sph[i].enuc += sph[i].dnuc;
    }
}

template <class Tsph>
void correctReleasedNuclearEnergy(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].cmps  = sph[i].cmps0;
        sph[i].dnuc  = GeneralSPH::calcReleasedNuclearEnergy(RP::Timestep*0.5,
                                                             sph[i].dens0,
                                                             sph[i].temp0,
                                                             sph[i].cmps.getPointer());
        sph[i].dnuc += GeneralSPH::calcReleasedNuclearEnergy(RP::Timestep*0.5,
                                                             sph[i].dens,
                                                             sph[i].temp,
                                                             sph[i].cmps.getPointer());
        sph[i].enuc += sph[i].dnuc;
    }
}
*/

template <class Tsph>
void calcBalsaraSwitch(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].calcBalsaraSwitch();
    }
}

template <class Tsph,
          class Tbhns>
void addAdditionalForce(Tsph & sph,
                        Tbhns & bhns) {
    if(RP::FlagDamping != 2) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].addAdditionalForce();
        }
    } else {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].addAdditionalForceDamping2();
        }
        for(PS::S64 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
            bhns[i].addAdditionalForceDamping2();
        }
    }
}

template <class Tsph>
void calcAlphaDot(Tsph & sph) {
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].calcAlphaDot();
    }
}

template <class Tsph>
void sumAcceleration(Tsph & sph) {
    if(RP::FlagGravity == 0) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].copyAcceleration();
        }
    } else {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].sumAcceleration();
        }
    }
}

template <class Tsph,
          class Tbhns>
PS::F64 calcPotentialEnergy(Tsph & sph,
                            Tbhns & bhns) {
    PS::F64 eloc = 0.;
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        eloc += sph[i].calcPotentialEnergy();
    }
    for(PS::S64 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
        eloc += bhns[i].calcPotentialEnergy();
    }
    PS::F64 eglb = PS::Comm::getSum(eloc);
    return eglb;
}

template <class Tsph,
          class Tbhns>
PS::F64 calcEnergy(Tsph & sph,
                   Tbhns & bhns) {
    PS::F64 eloc = 0.;
    if(RP::FlagDamping != 2) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            eloc += sph[i].calcEnergy();
        }
        for(PS::S64 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
            eloc += bhns[i].calcEnergy();
        }
    } else {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            eloc += sph[i].calcEnergyDamping2();
        }
        for(PS::S64 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
            eloc += bhns[i].calcEnergyDamping2();
        }
    }
    PS::F64 eglb = PS::Comm::getSum(eloc);
    return eglb;
}

template <class Tsph>
PS::F64 calcReleasedNuclearEnergyTotal(Tsph & sph) {
    PS::F64 eloc = 0.;
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        eloc += sph[i].mass * sph[i].enuc;
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
PS::F64 calcTimestep(Tsph & sph) {
    PS::S32 nloc = sph.getNumberOfParticleLocal();
    PS::F64 dtc = RP::MaximumTimestep;
    for(PS::S32 i = 0; i < nloc; i++) {
        PS::F64 dttmp = sph[i].calcTimestep();
        if(dttmp < dtc)
            dtc = dttmp;
    }
    dtc = PS::Comm::getMinValue(dtc);

    PS::F64 dt = RP::Timestep;
    if(dt > dtc) {
        while(dt > dtc) {
            dt *= 0.5;
            if(dt <= RP::MinimumTimestep) {
                if(PS::Comm::getRank() == 0) {
                    printf("Time: %.10f Timestep: %+e\n", RP::Time, dtc);
                }
                sph.writeParticleAscii("snap/smalldt.dat");
                PS::Comm::barrier();
                PS::Abort();
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

template <class Tsph>
void dumpHighEnergyParticle(Tsph & sph) {
    bool floc = false;
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        //PS::F64 umax = CalcEquationOfState::getEnergyMax(sph[i].dens, sph[i].abar, sph[i].zbar);
        PS::F64 umax = CalcEquationOfState::getEnergyMax(sph[i].dens, sph[i].cmps);
        if(sph[i].uene > umax) {
            floc = true;
        }
    }
    bool fglb = PS::Comm::synchronizeConditionalBranchOR(floc);
    if(fglb) {
        sph.writeParticleAscii("snap/hoge.dat");
        PS::Finalize();
        exit(0);
    }
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls,
          class Tdensity,
          class Thydro,
          class Tgravity>
void calcSPHKernel(Tdinfo & dinfo,
                   Tsph & sph,
                   Tbhns & bhns,
                   Tmsls & msls,
                   Tdensity & density,
                   Thydro & hydro,
                   Tgravity & gravity) {
//#define NBODYLIKE
#ifdef NBODYLIKE
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].acc   = 0.0;
        sph[i].acch  = 0.0;
        sph[i].accg1 = 0.0;
        sph[i].accg2 = 0.0;
        sph[i].pot   = 0.0;
        sph[i].pot3  = 0.0;        
    }
    addAdditionalForce(sph, bhns);
#else
    WT::start();
    RP::KernelSupportRadiusMaximum = calcKernelSupportRadiusMaximum(sph);
    calcDensityKernel(dinfo, sph, density);
    WT::accumulateCalcDensity();
    WT::start();
    calcAbarZbar(sph);
    referEquationOfState(sph);
    calcBalsaraSwitch(sph);
    WT::accumulateOthers();
    WT::start();
    calcGravityKernel(dinfo, sph, bhns, msls, gravity);
    WT::accumulateCalcGravity();
    WT::start();
    hydro.calcForceAllAndWriteBack(calcHydro(), sph, dinfo);
    WT::accumulateCalcHydro();
    WT::start();
    sumAcceleration(sph);
    addAdditionalForce(sph, bhns);
    calcAlphaDot(sph);
    WT::accumulateOthers();
#endif
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
    PS::ParticleSystem<MassLess> msls;
    msls.initialize();
    msls.createParticle(0);
    msls.setNumberOfParticleLocal(0);
    PS::ParticleSystem<BlackHoleNeutronStar> bhns;
    bhns.initialize();
    bhns.createParticle(0);
    bhns.setNumberOfParticleLocal(0);
    PS::TreeForForceShort<Density, DensityEPI, DensityEPJ>::Gather density;
    density.initialize(0);
    PS::TreeForForceShort<Hydro, HydroEPI, HydroEPJ>::Symmetry hydro;
    hydro.initialize(0);
    PS::TreeForForce<PS::SEARCH_MODE_LONG, Gravity, GravityEPI, GravityEPJ,
                     PS::GravityMonopole, PS::GravityMonopole, PS::GravitySPJ> gravity;
//  160610 from
    gravity.initialize(0, 0.5);
//    gravity.initialize(0, 0.4);
//  160610 to

    initializeSimulation();

    if(atoi(argv[1]) == 0) {
        startSimulation(argv, dinfo, sph, bhns, msls, density, hydro, gravity);
    } else {
        restartSimulation(argv, dinfo, sph, bhns, msls);
    }

    loopSimulation(dinfo, sph, bhns, msls, density, hydro, gravity);

    finalizeSimulation(dinfo, sph, bhns, msls);

    PS::Finalize();

    return 0;
}
