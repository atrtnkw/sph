#include <iostream>
#include <vector>
#include <cassert>

#include "particle_simulator.hpp"
#include "hdr_sph.hpp"

#include "hdr_kernel.hpp"
#include "hdr_density.hpp"
#include "hdr_hydro.hpp"

#include "gp5util.h"
#include "hdr_gravity.hpp"

static PS::S32 nptclmax = 65536;

template <class Theader>
void setParameterParticle(Theader & header);

template <class Tptcl>
void writeSnapshot(FILE *fp,
                   Tptcl & system);

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

    sph.adjustPositionIntoRootDomain(dinfo);
    dinfo.decomposeDomainAll(sph);

    density.initialize(nptclmax);
    derivative.initialize(nptclmax);

    sph.exchangeParticle(dinfo);
    calcSPHKernel(dinfo, sph, density, derivative,
                  calcDensity(), calcDerivative());

#ifdef GRAVITY
    PS::TreeForForceLong<Gravity, GravityEPI, GravityEPJ>::Monopole gravity;
    gravity.initialize(nptclmax);
    g5_open();
    g5_set_eps_to_all(SPH::eps);
    gravity.calcForceAllAndWriteBack(calcGravity<GravityEPJ>(),
                                     calcGravity<PS::SPJMonopole>(),
                                     sph,
                                     dinfo);
#endif

    dtime = calcTimeStep(sph, time, 1 / 64.);

    PS::F64 tout = 0.0;
    PS::F64 dtdc = 0.25;
    PS::S32 nstp = 0;
    while(time < tend){
        if(time >= tout) {
            char filename[64];
            sprintf(filename, "snap/sph_t%04d.dat", nstp);
            sph.writeParticleAscii(filename);
            tout += dtsp;
            nstp++;
        }
        
        PS::F64 etot = calcEnergy(sph);
        if(rank == 0) {
            fprintf(stderr, "time: %.10f %+e %+e\n", time, dtime, etot);
        }

        predict(sph, dtime);
        if(SPH::cbox.low_[0] != SPH::cbox.high_[0]) {
            sph.adjustPositionIntoRootDomain(dinfo);
        }

        if(time - (PS::S64)(time / dtdc) * dtdc == 0.){
            dinfo.decomposeDomainAll(sph);
        }

        sph.exchangeParticle(dinfo);
        calcSPHKernel(dinfo, sph, density, derivative,
                      calcDensity(), calcDerivative());

#ifdef GRAVITY
        gravity.calcForceAllAndWriteBack(calcGravity<GravityEPJ>(),
                                         calcGravity<PS::SPJMonopole>(),
                                         sph,
                                         dinfo);
#endif
        
        correct(sph, dtime);

        time += dtime;
        dtime = calcTimeStep(sph, time, dtime);
    }

    {
        char filename[64];
        sprintf(filename, "snap/sph_t%04d.dat", nstp);
        sph.writeParticleAscii(filename);
    }
    
    PS::Finalize();

    return 0;
}

template <class Theader>
void setParameterParticle(Theader & header) {
    SPH::cbox  = header.cbox;
    SPH::eta   = header.eta;
    SPH::cinv  = header.gamma - 1.d;
    SPH::alphamax = header.alphamax;
    SPH::alphamin = header.alphamin;
    SPH::tceff = header.tceff;
    
    return;
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
    referEquationOfState(sph);
    calcBalsaraSwitch(sph);
    derivative.calcForceAllAndWriteBack(calcDerivative, sph, dinfo);
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].calcAlphaDot();
    }

    return;
}

