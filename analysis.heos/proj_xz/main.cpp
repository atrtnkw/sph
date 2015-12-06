#include <iostream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstring>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};
class Quantity;

#include "particle_simulator.hpp"
#include "hdr_time.hpp"
#include "hdr_run.hpp"
#include "vector_x86.hpp"
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_agas.hpp"
#include "hdr_massless.hpp"
#include "hdr_quantity.hpp"

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = PS::Comm::getRank();
    PS::S32 size = PS::Comm::getNumberOfProc();
    Header hdr;
    PS::DomainInfo dinfo;
    dinfo.initialize();
    PS::ParticleSystem<GeneralSPH> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);
    PS::ParticleSystem<MassLess> msls;
    msls.initialize();
    msls.createParticle(0);
    PS::TreeForForceShort<Quantity, QuantityEPI, QuantityEPJ>::Scatter quantity;
    quantity.initialize(0);

    if(PS::Comm::getRank() == 0) {
        hdr.initialize(argv[1]);
    }
    PS::Comm::broadcast(&hdr, 1);
    sph.readParticleAscii(argv[2], hdr);

    convertCoordinate(hdr.cpos, sph);

    if(rank == 0) {
        generateMassLessParticlePlainXZ(hdr, msls);
    }

    dinfo.decomposeDomainAll(sph);
    sph.exchangeParticle(dinfo);
    msls.exchangeParticle(dinfo);

    quantity.setParticleLocalTree(sph);
    quantity.setParticleLocalTree(msls, false);
    quantity.calcForceMakingTree(calcQuantity(), dinfo);
    PS::S32 nsph  = sph.getNumberOfParticleLocal();
    PS::S32 nmsls = msls.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nmsls; i++) {
        msls[i].copyFromForce(quantity.getForce(i+nsph));
    }

    msls.writeParticleAscii(argv[3]);

    PS::Finalize();

    return 0;
}
