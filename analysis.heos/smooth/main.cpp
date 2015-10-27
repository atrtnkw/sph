#include <iostream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sys/time.h>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};

#include "particle_simulator.hpp"

#include "vector_x86.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_hgas.hpp"
#include "hdr_quantity.hpp"

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = PS::Comm::getRank();
    PS::S32 size = PS::Comm::getNumberOfProc();
    Header header;
    PS::DomainInfo dinfo;
    dinfo.initialize();
    PS::ParticleSystem<GeneralSPH> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);
    PS::TreeForForceShort<Quantity, QuantityEPI, QuantityEPJ>::Symmetry quantity;
    quantity.initialize(0);

    header.nptcl = atoi(argv[1]);
    sph.readParticleAscii(argv[2], header);

    dinfo.decomposeDomainAll(sph);
    sph.exchangeParticle(dinfo);

    quantity.calcForceAllAndWriteBack(calcQuantity(), sph, dinfo);

    sph.writeParticleAscii(argv[3]);

    PS::Finalize();

    return 0;
}
