#include <iostream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstring>

#include "particle_simulator.hpp"

static PS::F64 MinimumTimeStep    = 1. / 64.;
static PS::F64 MinimumTimeStepInv = 1. / MinimumTimeStep;

#include "hdr_time.hpp"
#include "vector_x86.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_massless.hpp"
#include "hdr_density.hpp"
#include "hdr_hydro.hpp"

#include "gp5util.h"
#include "hdr_gravity.hpp"

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

template <class Thdr,
          class Tdinfo,
          class Tpsys,
          class Ttree1,
          class Ttree2,
          class Tfunc1,
          class Tfunc2>
void calcSPHKernel(Thdr & header,
                   Tdinfo & dinfo,
                   Tpsys & sph,
                   Ttree1 & density,
                   Ttree2 & derivative,
                   Tfunc1 calcDensity,
                   Tfunc2 calcDerivative);

template <class Tdinfo,
          class Tpsys1,
          class Tpsys2,
          class Ttree>
void calcGravityKernel(Tdinfo & dinfo,
                       Tpsys1 & sph,
                       Tpsys2 & msls,
                       Ttree & gravity);

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = PS::Comm::getRank();
    PS::S32 size = PS::Comm::getNumberOfProc();
    PS::F64 dtime;
    Header header;
    PS::DomainInfo dinfo;
    dinfo.initialize();
    PS::ParticleSystem<SPH> sph;
    sph.initialize();
    sph.createParticle(nptclmax);
    sph.setNumberOfParticleLocal(0);
    PS::ParticleSystem<MassLess> msls;
    msls.initialize();
    msls.createParticle(nptclmax);
    PS::TreeForForceShort<Density, DensityEPI, DensityEPJ>::Gather            density;
    density.initialize(nptclmax);
    PS::TreeForForceShort<Derivative, DerivativeEPI, DerivativeEPJ>::Symmetry derivative;
    derivative.initialize(nptclmax);        
#ifdef GRAVITY
    PS::TreeForForceLong<Gravity, GravityEPI, GravityEPJ>::Monopole gravity;
    gravity.initialize(nptclmax, 0.5);
    g5_open();
#endif

    WT::clear();

    if(atoi(argv[1]) == 0) {
        
        sph.readParticleAscii(argv[2], header);
#ifdef WD_DAMPINGB
        generateMassLessParticle(msls, sph);
#endif
        PS::Comm::broadcast(&header, 1);
        setParameterParticle(header);
        
        if(SPH::cbox.low_[0] != SPH::cbox.high_[0]) {
            dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
            dinfo.setPosRootDomain(SPH::cbox.low_, SPH::cbox.high_);
        }
        
        WT::start();
        sph.adjustPositionIntoRootDomain(dinfo);
        msls.adjustPositionIntoRootDomain(dinfo);
        WT::accumulateOthers();
        WT::start();
        PS::MT::init_genrand(0);
        dinfo.decomposeDomainAll(sph);
        WT::accumulateDecomposeDomain();
        
        WT::start();
        calcFieldVariable(sph);
        WT::accumulateOthers();
        
        WT::start();
        sph.exchangeParticle(dinfo);
#ifdef WD_DAMPINGB
        msls.exchangeParticle(dinfo);
#endif
        WT::accumulateExchangeParticle();
        WT::start();
        SPH::ksrmax   = calcSystemSize(sph);
        header.ksrmax = SPH::ksrmax;
        WT::accumulateOthers();
        calcSPHKernel(header, dinfo, sph, density, derivative,
                      calcDensity(), calcDerivative());

#ifdef GRAVITY
        WT::start();
        g5_set_eps_to_all(SPH::eps);
        calcGravityKernel(dinfo, sph, msls, gravity);
        WT::accumulateCalcGravity();
#endif
    } else {
        char filename[64];
        sprintf(filename, "%s_p%06d.hexa", argv[2], rank);
        readRestartFile(filename, header, dinfo, sph);
        setParameterParticle(header);

        if(SPH::cbox.low_[0] != SPH::cbox.high_[0]) {
            dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
            dinfo.setPosRootDomain(SPH::cbox.low_, SPH::cbox.high_);
        }

#ifdef GRAVITY
        g5_set_eps_to_all(SPH::eps);
#endif
    }

    FILE *fplog = fopen("snap/time.log", "w");
    FILE *fptim = fopen("snap/prof.log", "w");
    PS::F64 tout = header.time;
    const PS::S32 nstp = 4;
    while(header.time < header.tend){
        /*
        if(header.time > 90.0341) {
            PS::TimeProfile wt_dens = density.getTimeProfile();
            char filename[64];
            sprintf(filename, "snap/wt_%04d.dat", rank);
            FILE *fp = fopen(filename, "w");
            fprintf(fp, "%4d", rank);
            fprintf(fp, " %+e", wt_dens.make_local_tree);
            fprintf(fp, " %+e", wt_dens.make_global_tree);
            fprintf(fp, " %+e", wt_dens.calc_force);
            fprintf(fp, " %+e", wt_dens.calc_moment_local_tree);
            fprintf(fp, " %+e", wt_dens.calc_moment_global_tree);
            fprintf(fp, " %+e", wt_dens.make_LET_1st);
            fprintf(fp, " %+e", wt_dens.make_LET_2nd);
            fprintf(fp, " %+e", wt_dens.exchange_LET_1st);
            fprintf(fp, " %+e", wt_dens.exchange_LET_2nd);
            fprintf(fp, "\n");
            fclose(fp);
        }
        density.clearTimeProfile();
        */

        doThisEveryTime(dtime, tout, header, dinfo, sph, msls, fplog, fptim);

        WT::start();
        dtime = calcTimeStep(sph, header.time, MinimumTimeStep);
        WT::accumulateOthers();

        WT::start();
        predict(sph, dtime);
        WT::accumulateIntegrateOrbit();
        WT::start();
        if(SPH::cbox.low_[0] != SPH::cbox.high_[0]) {
            sph.adjustPositionIntoRootDomain(dinfo);
        }
        WT::accumulateOthers();

        WT::start();
        if(header.istp == nstp) {
            PS::MT::init_genrand(0);
            dinfo.decomposeDomainAll(sph);
            header.istp = 0;
        }
        header.istp++;
        WT::accumulateDecomposeDomain();

        WT::start();
        sph.exchangeParticle(dinfo);
        WT::accumulateExchangeParticle();

        calcSPHKernel(header, dinfo, sph, density, derivative,
                      calcDensity(), calcDerivative());

#ifdef GRAVITY
        WT::start();
        calcGravityKernel(dinfo, sph, msls, gravity);
        WT::accumulateCalcGravity();
#endif

        WT::start();
        correct(sph, dtime);
        WT::accumulateIntegrateOrbit();

        header.time += dtime;
    }

    fclose(fplog);
    fclose(fptim);

    finalizeSimulation(header.time, header, dinfo, sph, msls);
    
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

template <class Thdr,
          class Tdinfo,
          class Tpsys,
          class Ttree1,
          class Ttree2,
          class Tfunc1,
          class Tfunc2>
void calcSPHKernel(Thdr & header,
                   Tdinfo & dinfo,
                   Tpsys & sph,
                   Ttree1 & density,
                   Ttree2 & derivative,
                   Tfunc1 calcDensity,
                   Tfunc2 calcDerivative)
{
    const PS::F64 expand  = 1.1;
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].rs = expand * sph[i].ksr;
    }

    PS::S32 cnt = 0;
    for(bool repeat = true; repeat == true;) {
        WT::start();
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
        WT::accumulateCalcDensity();
        repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
        cnt++;
    }
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

template <class Tdinfo,
          class Tpsys1,
          class Tpsys2,
          class Ttree>
void calcGravityKernel(Tdinfo & dinfo,
                       Tpsys1 & sph,
                       Tpsys2 & msls,
                       Ttree & gravity)
{
#ifdef WD_DAMPINGB
    gravity.setParticleLocalTree(sph);
    gravity.setParticleLocalTree(msls, false);
    gravity.calcForceMakingTree(calcGravity<GravityEPJ>(),
                                calcGravity<PS::SPJMonopole>(),
                                dinfo);
    PS::S32 nsph  = sph.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nsph; i++) {
        sph[i].copyFromForce(gravity.getForce(i));
    }
    PS::S32 nmsls = msls.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nmsls; i++) {
        msls[i].copyFromForce(gravity.getForce(i+nsph));
    }
#else
    gravity.calcForceAllAndWriteBack(calcGravity<GravityEPJ>(),
                                     calcGravity<PS::SPJMonopole>(),
                                     sph,
                                     dinfo);
#endif
}
