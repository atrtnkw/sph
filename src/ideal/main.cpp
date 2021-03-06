#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>

#include "particle_simulator.hpp"

#include "hdr_time.hpp"

#include "vector_x86.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_massless.hpp"
#include "hdr_density.hpp"
#include "hdr_hydro.hpp"
#include "hdr_iad.hpp"
#include "hdr_gradient.hpp"

#include "gp5util.h"
#include "hdr_gravity.hpp"
#include "hdr_vgravity.hpp"

#include "hdr_damp.hpp"

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
          class Ttree3,
          class Ttree4,
          class Tfunc1,
          class Tfunc2,
          class Tfunc3,
          class Tfunc4_1,
          class Tfunc4_2>
void calcSPHKernel(Tdinfo & dinfo,
                   Tpsys & sph,
                   Ttree1 & density,
                   Ttree2 & derivative,
                   Ttree3 & gradient,
                   Ttree4 & gravity,
                   Tfunc1 calcDensity,
                   Tfunc2 calcDerivative,
                   Tfunc3 calcGradient,
                   Tfunc4_1 calcGravityEPJ,
                   Tfunc4_2 calcGravitySPJ);

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
    PS::TreeForForceShort<Gradient, GradientEPI, GradientEPJ>::Gather         gradient;

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
    gradient.initialize(nptclmax);

    WT::start();
    calcFieldVariable(sph);
    WT::accumulateOthers();

    SPH::epsu = setEpsilonOfInternalEnergy(sph);

#ifdef SYMMETRIZED_GRAVITY
    PS::TreeForForce<PS::SEARCH_MODE_LONG, Gravity,
        SymmetrizedGravityEPI, SymmetrizedGravityEPJ, PS::MomentMonopoleSymmetrized,
        PS::MomentMonopoleSymmetrized, PS::SPJMonopoleSymmetrized> gravity;
    gravity.initialize(nptclmax);
#else
    PS::TreeForForceLong<Gravity, GravityEPI, GravityEPJ>::Monopole gravity;
    gravity.initialize(nptclmax);
    g5_open();
    g5_set_eps_to_all(SPH::eps);
#endif

    WT::start();
    sph.exchangeParticle(dinfo);
    WT::accumulateExchangeParticle();
    calcSPHKernel(dinfo, sph, density, derivative, gradient, gravity,
                  calcDensity(), calcDerivative(), calcGradient(),
                  calcGravityEPJ(), calcGravitySPJ());

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
        }

        PS::F64 etot = calcEnergy(sph);
        WT::reduceInterProcess();
        if(rank == 0) {
            //fprintf(stderr, "time: %.10f %+e %+e\n", time, dtime, etot);
            fprintf(fplog, "time: %.10f %+e %+e\n", time, dtime, etot);
            fflush(fplog);
            WT::dump(time, fptim);
            fflush(fptim);
        }
        WT::clear();

#if TIMETEST
        if(time > 0.d) {
            fprintf(stdout, "dens: %6d grdh: %6d hydr: %6d\n", ncalcdens, ncalcgrdh, ncalchydr);
            fprintf(stdout, "dens: %+e grdh: %+e hydr: %+e\n", tcalcdens, tcalcgrdh, tcalchydr);
            fprintf(stdout, "dens: %+e grdh: %+e hydr: %+e\n",
                    tcalcdens/ncalcdens, tcalcgrdh/ncalcgrdh, tcalchydr/ncalchydr);
            fprintf(stdout, "ninteract: %d %d\n", ninteract, nintrhydr);
            fprintf(stdout, "\n");
            PS::Finalize();
            exit(0);
        } else {
            tcalcdens = 0.d;
            tcalcgrdh = 0.d;
            tcalchydr = 0.d;
            ncalcdens = 0;
            ncalcgrdh = 0;
            ncalchydr = 0;
            ninteract = 0;
            nintrhydr = 0;
            //tcalcsqrt = 0.d;
            //tcalcothr = 0.d;
            //ncalcsqrt = 0;
            //ncalcothr = 0;
        }
#endif
        
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

        calcSPHKernel(dinfo, sph, density, derivative, gradient, gravity,
                      calcDensity(), calcDerivative(), calcGradient(),
                      calcGravityEPJ(), calcGravitySPJ());


        /*
        sph.writeParticleAscii("snap/hoge.dat");
        PS::Finalize();
        exit(0);
        */

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
    // check
    const PS::F64 dtmin = 1e-16;

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
            //assert(dt > dtmin);
            if(dt <= dtmin) {
                system.writeParticleAscii("snap/timestep_too_small.dat");
                PS::Abort();
            }
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
          class Ttree3,
          class Ttree4,
          class Tfunc1,
          class Tfunc2,
          class Tfunc3,
          class Tfunc4_1,
          class Tfunc4_2>
void calcSPHKernel(Tdinfo & dinfo,
                   Tpsys & sph,
                   Ttree1 & density,
                   Ttree2 & derivative,
                   Ttree3 & gradient,
                   Ttree4 & gravity,
                   Tfunc1   calcDensity,
                   Tfunc2   calcDerivative,
                   Tfunc3   calcGradient,
                   Tfunc4_1 calcGravityEPJ,
                   Tfunc4_2 calcGravitySPJ)
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
#ifdef GRAVITY
    WT::start();
    gravity.calcForceAllAndWriteBack(calcGravityEPJ,
                                     calcGravitySPJ,
                                     sph,
                                     dinfo);
    WT::accumulateCalcGravity();    
#endif
#ifdef INTEGRAL_APPROACH_DERIVATIVE
    WT::start();
    gradient.calcForceAllAndWriteBack(calcGradient, sph, dinfo);
    WT::accumulateCalcHydro();
#endif
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

