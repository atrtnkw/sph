#include <iostream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sys/time.h>

#include "particle_simulator.hpp"

#include "vector_x86.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_massless.hpp"
#include "hdr_quantity.hpp"

static PS::S32 nptclmax = 65536;

class WallclockTime {
    PS::F64 TimeStart;
    PS::F64 TimeFileInput;
    PS::F64 TimeFileOutput;
    PS::F64 TimeFDPS;
    PS::F64 TimeOther;

    WallclockTime() {}
    ~WallclockTime() {}
    WallclockTime(const WallclockTime & c);
    WallclockTime & operator = (const WallclockTime & c);
    static WallclockTime & getInstance() {
        static WallclockTime inst;
        return inst;
    }

    PS::F64 getWallclockTime() {
        struct timeval tv;
        gettimeofday(& tv, NULL);
        return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6);
    }
public:

    static void start() {
        getInstance().TimeStart = getInstance().getWallclockTime();
    }

    static void accumulateTimeFileInput() {
        WallclockTime & p = getInstance();
        p.TimeFileInput += p.getWallclockTime() - p.TimeStart;
    }
    
    static void accumulateTimeFileOutput() {
        WallclockTime & p = getInstance();
        p.TimeFileOutput += p.getWallclockTime() - p.TimeStart;
    }

    static void accumulateTimeFDPS() {
        WallclockTime & p = getInstance();
        p.TimeFDPS += p.getWallclockTime() - p.TimeStart;
    }

    static void accumulateTimeOther() {
        WallclockTime & p = getInstance();
        p.TimeOther += p.getWallclockTime() - p.TimeStart;
    }

    static void outputTimeProfile() {
        WallclockTime & p = getInstance();
        printf("TimeFileInput:  %e\n", p.TimeFileInput);
        printf("TimeFileOutput: %e\n", p.TimeFileOutput);
        printf("TimeFDPS:       %e\n", p.TimeFDPS);
        printf("TimeOther:      %e\n", p.TimeOther);
    }
};

typedef WallclockTime WT;

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

    WT::start();
    header.nptcl = atoi(argv[1]);
    sph.readParticleAscii(argv[2], header);
    WT::accumulateTimeFileInput();

    WT::start();
    if(rank == 0) {
        generateMassLessParticle(msls);    
    }
    WT::accumulateTimeOther();

    WT::start();
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
    WT::accumulateTimeFDPS();

    WT::start();
    msls.writeParticleAscii(argv[3]);
    WT::accumulateTimeFileOutput();

    WT::outputTimeProfile();

    PS::Finalize();

    return 0;
}
