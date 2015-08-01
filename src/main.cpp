#include <iostream>
#include <vector>
#include <cassert>

#include "particle_simulator.hpp"

#include "vector_x86.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_density.hpp"
#include "hdr_hydro.hpp"

#include "gp5util.h"
#include "hdr_gravity.hpp"

#include "hdr_damp.hpp"

#include "hdr_time.hpp"

static PS::S32 nptclmax = 65536;

template <class Tptcl>
void referEquationOfState(Tptcl & system);

template <class Tptcl>
void calcBalsaraSwitch(Tptcl & system);

template <class Tptcl>
PS::F64 calcTimeStep(Tptcl & system,
                     PS::F64 t,
                     PS::F64 dt);

template <class Tptcl>
PS::F64 calcEnergy(Tptcl & system);

template <class Tptcl>
void predict(Tptcl & system,
             PS::F64 dt);
             
template <class Tptcl>
void correct(Tptcl & system,
             PS::F64 dt);

template <class Tptcl>
void addAdditionalForce(Tptcl & system,
                        PS::F64 dt);

template <class Tdinfo,
          class Tpsys,
          class Ttree1,
          class Ttree2,
          class Tfunc1,
          class Tfunc2>
void calcSPHKernel(Tdinfo & dinfo,
                   Tpsys & sph,
                   Ttree1 & density,
                   Ttree2 & derivative,
                   Tfunc1 calcDensity,
                   Tfunc2 calcDerivative);

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = PS::Comm::getRank();
    PS::S32 size = PS::Comm::getNumberOfProc();
    PS::F64 time, dtime, tend, dtsp;
    Header header;
    PS::DomainInfo dinfo;
    PS::ParticleSystem<SPH> sph;
    PS::TreeForForceShort<Density, DensityEPI, DensityEPJ>::Gather            density;
    PS::TreeForForceShort<Derivative, DerivativeEPI, DerivativeEPJ>::Symmetry derivative;

    sph.initialize();
    sph.createParticle(nptclmax);
    sph.setNumberOfParticleLocal(0);

    sph.readParticleAscii(argv[1], header);
    PS::Comm::broadcast(&header, 1);
    setParameterParticle(header);
    time      = header.time;
    tend      = header.tend;
    dtsp      = header.dtsp;

    dinfo.initialize();
    if(SPH::cbox.low_[0] != SPH::cbox.high_[0]) {
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
        dinfo.setPosRootDomain(SPH::cbox.low_, SPH::cbox.high_);
    }

    WT::clear();

    WT::start();
    sph.adjustPositionIntoRootDomain(dinfo);
    WT::accumulateOthers();
    WT::start();
    dinfo.decomposeDomainAll(sph);
    WT::accumulateDecomposeDomain();

    density.initialize(nptclmax);
    derivative.initialize(nptclmax);

    WT::start();
    calcFieldVariable(sph);
    WT::accumulateOthers();

    WT::start();
    sph.exchangeParticle(dinfo);
    WT::accumulateExchangeParticle();
    calcSPHKernel(dinfo, sph, density, derivative,
                  calcDensity(), calcDerivative());

#ifdef GRAVITY
    PS::TreeForForceLong<Gravity, GravityEPI, GravityEPJ>::Monopole gravity;
    gravity.initialize(nptclmax);
    g5_open();
    g5_set_eps_to_all(SPH::eps);
    WT::start();
    gravity.calcForceAllAndWriteBack(calcGravity<GravityEPJ>(),
                                     calcGravity<PS::SPJMonopole>(),
                                     sph,
                                     dinfo);
    WT::accumulateCalcGravity();
#endif

    FILE *fplog = fopen("snap/time.log", "w");
    FILE *fptim = fopen("snap/prof.log", "w");
    PS::F64 tout = time;
    PS::F64 dtdc = 0.25;
    PS::S32 nstp = 0;
    while(time < tend){
        dtime = calcTimeStep(sph, time, 1 / 64.);

        if(time >= tout) {
            char filename[64];
            sprintf(filename, "snap/sph_t%04d.dat", nstp);
            sph.writeParticleAscii(filename);
            tout += dtsp;
            nstp++;
//            PS::Finalize();
//            exit(0);
        }
        
        PS::F64 etot = calcEnergy(sph);
        WT::reduceInterProcess();
        if(rank == 0) {
            fprintf(fplog, "time: %.10f %+e %+e\n", time, dtime, etot);
            fflush(fplog);
//            fprintf(stderr, "time: %.10f %+e %+e\n", time, dtime, etot);
            WT::dump(time, fptim);
            fflush(fptim);
        }
        WT::clear();

        WT::start();
        predict(sph, dtime);
        WT::accumulateIntegrateOrbit();
        WT::start();
        if(SPH::cbox.low_[0] != SPH::cbox.high_[0]) {
            sph.adjustPositionIntoRootDomain(dinfo);
        }
        WT::accumulateOthers();

        WT::start();
        if(time - (PS::S64)(time / dtdc) * dtdc == 0.){
            dinfo.decomposeDomainAll(sph);
        }
        WT::accumulateDecomposeDomain();

        WT::start();
        sph.exchangeParticle(dinfo);
        WT::accumulateExchangeParticle();

        WT::start();
        calcFieldVariable(sph);
        WT::accumulateOthers();

        calcSPHKernel(dinfo, sph, density, derivative,
                      calcDensity(), calcDerivative());

        WT::start();
#ifdef GRAVITY
        gravity.calcForceAllAndWriteBack(calcGravity<GravityEPJ>(),
                                         calcGravity<PS::SPJMonopole>(),
                                         sph,
                                         dinfo);
#endif
        WT::accumulateCalcGravity();

        WT::start();
        correct(sph, dtime);
        WT::accumulateIntegrateOrbit();

#ifdef DAMPING
        DampVelocity::dampVelocity(sph, dtime);
        if(DampVelocity::stopDamping(time, sph)) {
            break;
        }
#endif

        time += dtime;
    }

    fclose(fplog);
    fclose(fptim);

    finalizeSimulation(nstp, sph);
    
    PS::Finalize();

    return 0;
}

template <class Tptcl>
void referEquationOfState(Tptcl & system) {
    
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++)
        system[i].referEquationOfState();

    return;
}

template <class Tptcl>
void calcBalsaraSwitch(Tptcl & system) {    
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++)
        system[i].calcBalsaraSwitch();

    return;
}

template <class Tptcl>
PS::F64 calcTimeStep(Tptcl & system,
                     PS::F64 t,
                     PS::F64 dt) {
    const PS::F64 dtmin = 1e-10;

    PS::S32 nloc = system.getNumberOfParticleLocal();
    PS::F64 dtc = 1e30;
    for(PS::S32 i = 0; i < nloc; i++) {
        PS::F64 dttmp = system[i].calcTimeStep();
        if(dttmp < dtc)
            dtc = dttmp;
    }
    dtc = PS::Comm::getMinValue(dtc);

    if(dt > dtc) {
        while(dt > dtc) {
            dt *= 0.5;
            assert(dt > dtmin);
        }        
    } else if(2. * dt <= dtc) {
        PS::F64 dt2 = 2. * dt;
        if(t - (PS::S64)(t / dt2) * dt2 == 0.) {
            dt *= 2.;
        }
    }

    return dt;
}

template <class Tptcl>
PS::F64 calcEnergy(Tptcl & system) {
    PS::F64 etot;
    PS::F64 eloc = 0.;
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0 ; i < nloc; i++)
        eloc += system[i].calcEnergy();
    etot = PS::Comm::getSum(eloc);

    return etot;
}

template <class Tptcl>
void predict(Tptcl & system,
             PS::F64 dt) {
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0 ; i < nloc; i++)
        system[i].predict(dt);    

    return;
}

template <class Tptcl>
void correct(Tptcl & system,
             PS::F64 dt) {
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0 ; i < nloc; i++)
        system[i].correct(dt);    

    return;
}

template <class Tptcl>
void addAdditionalForce(Tptcl & system) {
    PS::S32 nloc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0 ; i < nloc; i++)
        system[i].addAdditionalForce();    
    
    return;
}

template <class Tdinfo,
          class Tpsys,
          class Ttree1,
          class Ttree2,
          class Tfunc1,
          class Tfunc2>
void calcSPHKernel(Tdinfo & dinfo,
                   Tpsys & sph,
                   Ttree1 & density,
                   Ttree2 & derivative,
                   Tfunc1 calcDensity,
                   Tfunc2 calcDerivative)
{
    const PS::F64 expand = 1.1;
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].rs = expand * sph[i].ksr;
    }

    WT::start();
    PS::S32 cnt = 0;
    for(bool repeat = true; repeat == true;) {
        bool repeat_loc = false;
        repeat = false;
        density.calcForceAll(calcDensity, sph, dinfo);    
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].rs != 0.0) {
                if(density.getForce(i).itr == true) {
                    repeat_loc = true;
                    sph[i].rs *= expand;
                } else {
                    sph[i].rs = 0.0;
                    sph[i].copyFromForce(density.getForce(i));
                }
            }
        }
        repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
        cnt++;

    }
    WT::accumulateCalcDensity();
    WT::start();
    referEquationOfState(sph);
    WT::accumulateReferEquationOfState();
    WT::start();
    calcBalsaraSwitch(sph);
    WT::accumulateOthers();
    WT::start();
    derivative.calcForceAllAndWriteBack(calcDerivative, sph, dinfo);
    WT::accumulateCalcHydro();
    WT::start();
    addAdditionalForce(sph);
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].calcAlphaDot();
    }
    WT::accumulateOthers();

    return;
}

