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
#include "hdr_volume.hpp"
#include "hdr_auxiliary.hpp"
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

template <class Tdinfo,
          class Tsph,
          class Tvolume,
          class Tauxiliary>
void calcSPHKernel(Tdinfo & dinfo,
                   Tsph & sph,
                   Tvolume & volume,
                   Tauxiliary &auxiliary) {
    calcVolumeKernel(dinfo, sph, volume);
    auxiliary.calcForceAllAndWriteBack(calcAuxiliary(), sph, dinfo);
    calcBalsaraSwitch(sph);
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
    PS::TreeForForceShort<Volume, VolumeEPI, VolumeEPJ>::Gather volume;
    volume.initialize(0);
    PS::TreeForForceShort<Auxiliary, AuxiliaryEPI, AuxiliaryEPJ>::Gather auxiliary;
    auxiliary.initialize(0);

    if(atoi(argv[1]) == 0) {
        startSimulation(argv, dinfo, sph, volume, auxiliary);
    } else {
        restartSimulation(argv, dinfo, sph);
    }
    sph.writeParticleAscii("snap/hoge.dat");

    PS::Finalize();

    return 0;
}
