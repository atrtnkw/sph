#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>
#include <mpi.h>

enum KernelType {CubicSpline = 0, WendlandC2 = 1, WendlandC4 = 2};

#include "particle_simulator.hpp"
#include "hdr_time.hpp"
#include "hdr_run.hpp"
#include "vector_x86.hpp"
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_hgas.hpp"

extern "C" {
    void coldwd_(PS::F64 *rhoc, PS::S32 *ncolumn, PS::F64 *radius,
                 PS::F64 *density, PS::F64 *massrad, PS::F64 *totalmass);
}

PS::F64 makeDensityDistributionOfWhiteDwarf(PS::F64 totalmass,
                                            PS::S32 & ncolumn,
                                            PS::F64 * radius,
                                            PS::F64 * density,
                                            PS::F64 * massrad) {
    const PS::F64 tolerance = 0.0001;
    PS::F64 totalmasstemp = 0.;
    PS::F64 rhocmin = 1e4;
    PS::F64 rhocmax = 1e10;
    PS::S32 ncount  = 0;

    do {
        assert(ncount < 1000);
        ncount++;
        PS::F64 rhoc = sqrt(rhocmin * rhocmax);
        coldwd_(&rhoc, &ncolumn, radius, density, massrad, &totalmasstemp);
        if(totalmasstemp < totalmass) {
            rhocmin = rhoc;
        } else {
            rhocmax = rhoc;
        }
    }while(fabs(totalmasstemp - totalmass) > totalmass * tolerance);

    fprintf(stderr, "Totalmass: %+e, rhoc %+e, ncount: %8d\n", totalmasstemp, rhocmax, ncount);
    
    return totalmasstemp * CodeUnit::SolarMass;
}

template <class Tsph>
void makeHexagonalClosePackingOfUnitSphere(PS::S64 nptcl,
                                           PS::F64 totalmass,
                                           Tsph & sph) {
    PS::S64 nsearch  = 2 * nptcl;
    PS::F64 radius   = cbrt(sqrt(2.)*M_PI/(6.*nptcl));
    PS::F64 diameter = 2. * radius;
    PS::S32 nx = (PS::S32)cbrt(nsearch) + 1;
    PS::S32 ny = (PS::S32)cbrt(nsearch) + 1;
    PS::S32 nz = (PS::S32)cbrt(nsearch) + 1;
    PS::F64 dx = diameter;
    PS::F64 dy = sqrt(3.)*0.5*diameter;
    PS::F64 dz = sqrt(3.)*0.5*diameter;

    std::vector<PS::F64vec> ppos;
    PS::S32 npos = 0;

    for(PS::S32 iz = 0; iz < nz; iz++) {
        PS::F64vec tpos;
        tpos[2] = -1. + dz * iz;
        PS::F64 pz = 0.;
        for(PS::S32 iy = 0; iy < ny; iy++) {
            tpos[1] = -1. + dy * iy;
            for(PS::S32 ix = 0; ix < nx; ix++) {
                tpos[0] = -1. + dx * ix + (0.5 * dx) * (iy % 2) + (0.5 * dx) * (iz % 2);
                if(tpos * tpos < 1.) {
                    ppos.push_back(tpos);
                    npos++;
                }
            }
        }
    }

    PS::F64 mass = totalmass / (PS::F64)npos;
    sph.setNumberOfParticleLocal(npos);
    for(PS::S32 i = 0; i < npos; i++) {
        sph[i].id   = i;
        sph[i].mass = mass;
        sph[i].pos  = ppos[i];
    }

    fprintf(stderr, "TotalNumberOfParticle: %8lld\n", nptcl);

    return;
}

template <class Tsph>
void extendUnitSphere(PS::F64 totalmass,
                      PS::S32 ncolumn,
                      PS::F64 * radius,
                      PS::F64 * density,
                      PS::F64 * massrad,
                      Tsph & sph) {
    PS::F64 radius0[ncolumn];
    for(PS::S32 i = 0; i < ncolumn; i++) {
        radius0[i] = cbrt(massrad[i] / totalmass);
//        printf("hoge %+e %+e %+e %+e\n", radius0[i], radius[i], massrad[i], totalmass);
    }

    PS::F64 radfornr0[ncolumn+1];
    for(PS::S32 i = 0; i < ncolumn+1; i++) {
        if(i == 0) {
            radfornr0[i] = 0.;
        } else {
            radfornr0[i] = radius0[i+1];
        }
    }

    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        PS::F64 rptcl0 = sqrt(sph[i].pos * sph[i].pos);
        PS::S64 imin0  = 0;
        PS::S64 imax0  = ncolumn;
        PS::S64 imid0  = 0;
        PS::S64 ncount = 0;
        do {
            imid0 = (imin0 + imax0) / 2;
            if(rptcl0 < radfornr0[imid0]) {
                imax0 = imid0;
            } else {
                imin0 = imid0;
            }
            ncount++;
            assert(ncount < 1000);
            assert(imax0 != imin0);
        } while(imax0 - imin0 > 1);
        assert(imax0 > 0);
        PS::F64 facdr0 = (radfornr0[imax0] - rptcl0) / (radfornr0[imax0] - radfornr0[imax0-1]);
        PS::F64 rmax   = radius[imax0-1];
        PS::F64 rmin   = ((imax0 == 1) ? 0. : radius[imax0-2]);
        PS::F64 rptcl  = rmax - (rmax - rmin) * facdr0;
        PS::F64vec tpos = sph[i].pos * (rptcl / rptcl0);
        sph[i].id    = i;
        sph[i].istar = 0;
        sph[i].mass  = sph[i].mass * CodeUnit::UnitOfMassInv;
        sph[i].pos   = tpos * CodeUnit::UnitOfLengthInv;
        sph[i].vel   = 0.;
        sph[i].alph  = 1.;
        sph[i].alphu = 0.;
        sph[i].cmps[1] = 0.5;
        sph[i].cmps[2] = 0.5;
        sph[i].dens  = density[imax0-1] * CodeUnit::UnitOfDensityInv;
        sph[i].ksr   = SK::eta * SK::ksrh * cbrt(sph[i].mass / sph[i].dens);
        sph[i].uene  = CalcEquationOfState::getEnergyGivenTemperature(sph[i].dens, 1e7,
                                                                      sph[i].cmps);
    }

    return;
}

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    PS::S64 nptcl;
    PS::F64 totalmasspermsun;

    FILE *fp = fopen("input.list", "r");
    fscanf(fp, "%ld%lf", &nptcl, &totalmasspermsun);
    fclose(fp);
    fprintf(stderr, "N: %ld M: %.3f\n", nptcl, totalmasspermsun);

    const PS::S32 ncolumnmax = 600;
    PS::S32 ncolumn    = 0;
    PS::F64 radius[ncolumnmax], density[ncolumnmax], massrad[ncolumnmax];
    PS::F64 totalmass = makeDensityDistributionOfWhiteDwarf(totalmasspermsun, ncolumn, radius,
                                                            density, massrad);
    PS::ParticleSystem<GeneralSPH> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);
    makeHexagonalClosePackingOfUnitSphere(nptcl, totalmass, sph);

    init_flash_helmholtz_(&CodeUnit::FractionOfCoulombCorrection);
    extendUnitSphere(totalmass, ncolumn, radius, density, massrad, sph);

    sph.writeParticleAscii("fuga.dat");


    MPI_Finalize();

    return 0;
}
