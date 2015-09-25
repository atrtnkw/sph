#include <iostream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstring>

#include "particle_simulator.hpp"

#include "vector_x86.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_massless.hpp"
#include "hdr_quantity.hpp"

static PS::S32 nptclmax = 65536;

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = PS::Comm::getRank();
    PS::S32 size = PS::Comm::getNumberOfProc();
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
    PS::TreeForForceShort<Quantity, QuantityEPI, QuantityEPJ>::Scatter quantity;
    quantity.initialize(nptclmax);

    sph.readParticleAscii(argv[1], header);
    if(rank == 0) {
        generateMassLessParticle(msls);    
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

    msls.writeParticleAscii(argv[2]);

    PS::Finalize();

    return 0;
}
