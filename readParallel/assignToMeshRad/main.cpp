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
#ifdef USE_INTRINSICS
#include "vector_x86.hpp"
#endif
#include "hdr_dimension.hpp"
#include "hdr_kernel.hpp"
#include "hdr_sph.hpp"
#include "hdr_hgas.hpp"
#include "hdr_bhns.hpp"

class SPHDetailElement : public HelmholtzGas {
public:

    static       PS::S64 CompanionCO;
    static       PS::F64 VelocityOverRadius;
    static const PS::S64 NumberOfElement = 34;
    PS::F64 elem[NumberOfElement];

    SPHDetailElement() {
        this->id    = 0;
        this->istar = 0;
        this->mass  = 0.;
        this->pos   = 0.;
        this->vel   = 0.;
        this->dens  = 0.;
        this->ksr   = 0.;
        this->pot   = 0.;
        for(PS::S64 k = 0; k < NumberOfElement; k++) {
            this->elem[k] = 0.;
        }
    }

    void readAscii(FILE * fp) {
        fscanf(fp, "%lld%lld%lf", &this->id,     &this->istar,  &this->mass);   //   3
        fscanf(fp, "%lf%lf%lf",   &this->pos[0], &this->pos[1], &this->pos[2]); //   6
        fscanf(fp, "%lf%lf%lf",   &this->vel[0], &this->vel[1], &this->vel[2]); //   9
        fscanf(fp, "%lf%lf%lf",   &this->dens,   &this->ksr,    &this->pot);    //  12
        if(CompanionCO == 1) {
            this->elem[9]  = 0.49179301803663630124; // C
            this->elem[11] = 0.49179301803663630124; // O
            this->elem[13] = 0.01342153152708231446; // Ne
            this->elem[15] = 0.00070522296620528707; // Mg
            this->elem[17] = 0.00066876453890332400; // Si
            this->elem[19] = 0.00031136169559311488; // S
            this->elem[29] = 0.00130708319894335710; // Fe
        } else if(CompanionCO == 0) {
            this->elem[5]  = 0.59232608581213036673; // He
            this->elem[9]  = 0.19671720721465452049; // C
            this->elem[10] = 0.00512458476488597461; // N
            this->elem[11] = 0.19671720721465452049; // O
            this->elem[13] = 0.00075386998319660882 + 0.00536861261083292578; // Ne
            this->elem[15] = 0.00070522296620528707; // Mg
            this->elem[17] = 0.00066876453890332400; // Si
            this->elem[19] = 0.00031136169559311488; // S
            this->elem[29] = 0.00130708319894335710; // Fe
        } else {
            assert(NULL);
        }
        if(this->istar == 0) {
            PS::S64 dum;
            fscanf(fp, "%lld", &dum);
            for(PS::S64 k = 0; k < NumberOfElement; k++) {
                fscanf(fp, "%lf", &this->elem[k]);
            }
        }
    }

};
PS::S64 SPHDetailElement::CompanionCO        = -1;
PS::F64 SPHDetailElement::VelocityOverRadius = -1.;

namespace AssignToMesh {
#if 0
    const PS::S64 NumberOfRadius      = 128;
    const PS::S64 NumberOfInclination = 64;
    const PS::S64 NumberOfAzimuth     = 64;
#else
    const PS::S64 NumberOfRadius      = 64;
    const PS::S64 NumberOfInclination = 16;
    const PS::S64 NumberOfAzimuth     = 16;
#endif
    const PS::S64 NumberOfMesh3       = NumberOfRadius
        * NumberOfInclination * NumberOfAzimuth;
    const PS::F64 RadiusOfBox    = 2.4e11;
    const PS::F64 DeltaRadius    = RadiusOfBox / ((PS::F64)NumberOfRadius);
    const PS::F64 DeltaRadiusInv = 1. / DeltaRadius;
    static PS::F32 dens_l[NumberOfAzimuth][NumberOfInclination][NumberOfRadius];
    static PS::F32 dens_g[NumberOfAzimuth][NumberOfInclination][NumberOfRadius];
    static PS::F32 temp_l[NumberOfAzimuth][NumberOfInclination][NumberOfRadius];
    static PS::F32 temp_g[NumberOfAzimuth][NumberOfInclination][NumberOfRadius];
    static PS::F32 cmps_l[SPHDetailElement::NumberOfElement][NumberOfAzimuth][NumberOfInclination][NumberOfRadius];
    static PS::F32 cmps_g[SPHDetailElement::NumberOfElement][NumberOfAzimuth][NumberOfInclination][NumberOfRadius];
    static bool FirstCallInclination = true;
    static PS::F64 incl[NumberOfInclination];

    PS::S64 getIndexRadius(PS::F64 r) {
        PS::F64 fx = r * DeltaRadiusInv + 0.5;
        PS::S64 ix = (PS::S64)fx;
        ix = (ix >= 0)             ? ix : 0;
        ix = (ix < NumberOfRadius) ? ix : (NumberOfRadius - 1);
        return ix;
    }

    PS::F64 getPositionRadius(PS::S64 index) {
        return (DeltaRadius * ((PS::F64)index + 0.5));
    }

    void calcInclinationGrid() {
        if(FirstCallInclination) {
            FirstCallInclination = false;
            incl[0] = acos(1. - 2. / (PS::F64)NumberOfInclination);
            for(PS::S64 k = 1; k < NumberOfInclination - 1; k++) {
                incl[k] = acos(cos(incl[k-1]) - 2. / (PS::F64)NumberOfInclination);
            }
            incl[NumberOfInclination-1] = M_PI;
        }
    }

    PS::S64 getIndexInclination(PS::F64 theta) {
        calcInclinationGrid();
        PS::S64 index = -1;
        if(theta <= incl[0]) {
            index = 0;
        } else if (theta >= incl[NumberOfInclination-1]) {
            index = NumberOfInclination-1;
        } else {
                PS::S64 idn = 0;
                PS::S64 iup = NumberOfInclination - 1;
            do {
                PS::S64 imd = (idn + iup) / 2;
                if(theta < incl[imd]) {
                    iup = imd;
                } else {
                    idn = imd;
                }
                //assert(idn <= iup);
                assert(idn < iup);
            } while(iup - idn > 1);
            index = iup;
        }
        assert(0 <= index && index < NumberOfInclination);
        return index;
    }

    PS::F64 getPositionInclination(PS::S64 index) {
        calcInclinationGrid();
        return incl[index];
    }

    PS::F64 getPositionAzimuth(PS::S64 index) {
        PS::F64 delta = 2. * M_PI / (PS::F64)NumberOfAzimuth;
        return (delta * (PS::F64)(index+1));
    }

}

class SPHAnalysis : public HelmholtzGas {
public:
    SPHAnalysis() {
        this->id    = 0;
        this->istar = 0;
        this->mass  = 0.;
        this->pos   = 0.;
        this->vel   = 0.;
        this->uene  = 0.;
        this->alph  = 0.;
        this->alphu = 0.;
        this->ksr   = 0.;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            this->cmps[k] = 0.;
        }        
    }

    void readAscii(FILE * fp) {
        fscanf(fp, "%lld%lld%lf", &this->id, &this->istar, &this->mass);      //  3
        fscanf(fp, "%lf%lf%lf", &this->pos[0], &this->pos[1], &this->pos[2]); //  6
        fscanf(fp, "%lf%lf%lf", &this->vel[0], &this->vel[1], &this->vel[2]); //  9
        fscanf(fp, "%lf%lf%lf", &this->acc[0], &this->acc[1], &this->acc[2]); // 12
        fscanf(fp, "%lf%lf%lf", &this->uene, &this->alph, &this->alphu);      // 15
        fscanf(fp, "%lf%lf%lld", &this->dens, &this->ksr,  &this->np);        // 18
        fscanf(fp, "%lf%lf%lf", &this->vsnd, &this->pres, &this->temp);       // 21
        fscanf(fp, "%lf%lf%lf", &this->divv, &this->rotv, &this->bswt);       // 24
        fscanf(fp, "%lf%lf%lf", &this->pot,  &this->abar, &this->zbar);       // 27
        fscanf(fp, "%lf",       &this->enuc);                                 // 28
        fscanf(fp, "%lf%lf%lf", &this->vsmx, &this->udot, &this->dnuc);       // 31
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) {       // 32 -- 44
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        fscanf(fp, "%lf", &this->pot3);
        fscanf(fp, "%lf%lf%lf", &this->tempmax[0], &this->tempmax[1], &this->tempmax[2]);
        fscanf(fp, "%lf", &this->entr);
    }

};

template <class Tsph>
void obtainDensityCenter(Tsph & sph,
                         PS::F64vec & pcenter,
                         PS::F64vec & vcenter) {
    PS::S64    np   = 0;
    PS::F64    dloc = 0.;
    PS::F64vec ploc(0.);
    PS::F64vec vloc(0.);
    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        if(sph[ip].istar == 0) {
            continue;
        }
        np++;
        dloc += sph[ip].dens;
        ploc += sph[ip].dens * sph[ip].pos;
        vloc += sph[ip].dens * sph[ip].vel;
    }

    PS::F64    dglb;
    PS::F64vec pglb;
    PS::F64vec vglb;
    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(&dloc, &dglb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&ploc, &pglb, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&vloc, &vglb, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    pcenter = pglb / dglb;
    vcenter = vglb / dglb;
}

template <class Tsph,
          class Tsphelem>
void assignToMesh(char * ofile,
                  Tsph & sph,
                  Tsphelem & sphelem) {
    using namespace AssignToMesh;

    PS::S64 NumberOfElement = SPHDetailElement::NumberOfElement;

    PS::F64vec pcenter(0.);
    PS::F64vec vcenter(0.);
    obtainDensityCenter(sph, pcenter, vcenter);

    for(PS::S64 iz = 0; iz < NumberOfAzimuth; iz++) {
        for(PS::S64 iy = 0; iy < NumberOfInclination; iy++) {
            for(PS::S64 ix = 0; ix < NumberOfRadius; ix++) {
                dens_l[iz][iy][ix] = 0.;
                temp_l[iz][iy][ix] = 0.;
                for(PS::S64 k = 0; k < NumberOfElement; k++) {
                    cmps_l[k][iz][iy][ix] = 0.;
                }
            }
        }
    }

    for(PS::S64 ip = 0; ip < sph.getNumberOfParticleLocal(); ip++) {
        PS::F64vec pptcl = sph[ip].pos - pcenter;
        PS::F64vec vptcl = sph[ip].vel - vcenter;
        PS::F64    eptcl = 0.5 * (vptcl * vptcl) + sph[ip].pot;
        if(eptcl < 0.) {
            continue;
        }

        PS::F64 radius = sqrt(pptcl * pptcl);
        PS::F64 theta  = acos(pptcl[2] / radius);
        PS::F64 dtheta = (sph[ip].ksr <= radius) ? asin(sph[ip].ksr / radius) : M_PI;
        PS::S64 irmin  = getIndexRadius(radius-sph[ip].ksr);
        PS::S64 irmax  = getIndexRadius(radius+sph[ip].ksr);
        PS::S64 itmin  = getIndexInclination(theta-dtheta);
        PS::S64 itmax  = getIndexInclination(theta+dtheta);
        for(PS::S64 iphi = 0; iphi < NumberOfAzimuth; iphi++) {
            for(PS::S64 iinc = itmin; iinc <= itmax; iinc++) {
                for(PS::S64 irad = irmin; irad <= irmax; irad++) {
                    PS::F64 rad = getPositionRadius(irad);
                    PS::F64 inc = getPositionInclination(iinc);
                    PS::F64 phi = getPositionAzimuth(iphi);
                    PS::F64vec pos(rad * sin(inc) * cos(phi),
                                   rad * sin(inc) * sin(phi),
                                   rad * cos(inc));
                    PS::F64vec dr = pos - pptcl;
                    PS::F64 dr2   = dr * dr;
                    PS::F64 dr1   = sqrt(dr2);
                    PS::F64 hinv1 = 1. / sph[ip].ksr;
                    PS::F64 hinv3 = ND::calcVolumeInverse(hinv1);
                    PS::F64 dq1   = dr1 * hinv1;
                    PS::F64 kw0   = hinv3 * SK::kernel0th(dq1);
                    dens_l[iphi][iinc][irad] += sph[ip].mass * kw0;
                    PS::F64 dinv = 1. / sph[ip].dens;
                    PS::F64 ksph = sph[ip].mass * dinv * kw0;
                    temp_l[iphi][iinc][irad] += ksph * sph[ip].temp;
                }            
            }            
        }
    }

    for(PS::S64 ip = 0; ip < sphelem.getNumberOfParticleLocal(); ip++) {
        PS::F64vec pptcl = sphelem[ip].pos - pcenter;
        PS::F64vec vptcl = sphelem[ip].vel - vcenter;
        PS::F64    eptcl = 0.5 * (vptcl * vptcl) + sphelem[ip].pot;
        if(eptcl < 0.) {
            continue;
        }
        if(sphelem[ip].id % 320 != 0) {
            continue;
        }

        PS::F64 fexp  = pow(320., 1./3.);
        PS::F64 radius = sqrt(pptcl * pptcl);
        PS::F64 theta  = acos(pptcl[2] / radius);
        PS::F64 dtheta = (sphelem[ip].ksr*fexp <= radius) ?
            asin(sphelem[ip].ksr*fexp / radius) : M_PI;
        PS::S64 irmin  = getIndexRadius(radius-sphelem[ip].ksr*fexp);
        PS::S64 irmax  = getIndexRadius(radius+sphelem[ip].ksr*fexp);
        PS::S64 itmin  = getIndexInclination(theta-dtheta);
        PS::S64 itmax  = getIndexInclination(theta+dtheta);

        for(PS::S64 iphi = 0; iphi < NumberOfAzimuth; iphi++) {
            for(PS::S64 iinc = itmin; iinc <= itmax; iinc++) {
                for(PS::S64 irad = irmin; irad <= irmax; irad++) {
                    PS::F64 rad = getPositionRadius(irad);
                    PS::F64 inc = getPositionInclination(iinc);
                    PS::F64 phi = getPositionAzimuth(iphi);
                    PS::F64vec pos(rad * sin(inc) * cos(phi),
                                   rad * sin(inc) * sin(phi),
                                   rad * cos(inc));
                    PS::F64vec dr = pos - pptcl;
                    PS::F64 dr2   = dr * dr;
                    PS::F64 dr1   = sqrt(dr2);
                    PS::F64 hinv1 = 1. / (sphelem[ip].ksr * fexp);
                    PS::F64 dq1   = dr1 * hinv1;
                    PS::F64 hinv3 = ND::calcVolumeInverse(hinv1);
                    PS::F64 kw0   = hinv3 * SK::kernel0th(dq1);
                    PS::F64 dinv = 1. / sphelem[ip].dens;
                    PS::F64 ksph = sphelem[ip].mass * dinv * kw0;
                    for(PS::S64 k = 0; k < NumberOfElement; k++) {
                        cmps_l[k][iphi][iinc][irad] += ksph * sphelem[ip].elem[k];
                    }
                }            
            }            
        }
    }

    PS::S64 ierr = 0;
    ierr = MPI_Allreduce(dens_l, dens_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(temp_l, temp_g, NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    for(PS::S64 k = 0; k < NumberOfElement; k++) {
        ierr = MPI_Allreduce(cmps_l[k], cmps_g[k], NumberOfMesh3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }

    if(PS::Comm::getRank() == 0) {
#define MaedaFormat
#ifdef MaedaFormat
        using namespace AssignToMesh;
        FILE * fp = fopen(ofile, "w");
        fprintf(fp, "%d\n", 2);   // just flag
        fprintf(fp, "%+e\n", 1.); // time in day
        fprintf(fp, "%+e\n", 1.); // the same as above
        fprintf(fp, "%lld %lld %lld\n", NumberOfRadius, NumberOfInclination, NumberOfAzimuth);
        fprintf(fp, "\n");
        for(PS::S64 irad = 0; irad < NumberOfRadius; irad++) {
            fprintf(fp, "%+e\n", getPositionRadius(irad) / SPHDetailElement::VelocityOverRadius * 1e-5); // [km/s]
        }
        fprintf(fp, "\n");
        for(PS::S64 iinc = 0; iinc < NumberOfInclination; iinc++) {
            fprintf(fp, "%+e\n", getPositionInclination(iinc));
        }
        fprintf(fp, "\n");
        for(PS::S64 iphi = 0; iphi < NumberOfAzimuth; iphi++) {
            fprintf(fp, "%+e\n", getPositionAzimuth(iphi));
        }
        fprintf(fp, "\n");
        for(PS::S64 iphi = 0; iphi < NumberOfAzimuth; iphi++) {
            for(PS::S64 iinc = 0; iinc < NumberOfInclination; iinc++) {
                for(PS::S64 irad = 0; irad < NumberOfRadius; irad++) {
                    PS::F64 dens = (dens_g[iphi][iinc][irad] > 0.) ?
                        log10(dens_g[iphi][iinc][irad]) : -100.;
                    PS::F32 norm = 0.;
                    for(PS::S64 k = 0; k < NumberOfElement; k++) {
                        norm += cmps_g[k][iphi][iinc][irad];
                    }
                    PS::F32 ninv = (norm != 0. && dens_g[iphi][iinc][irad] != 0.) ? (1. / norm) : 0.;

                    fprintf(fp, "%4lld %4lld %4lld", irad+1, iinc+1, iphi+1);
                    fprintf(fp, " %+e %+e", dens, temp_g[iphi][iinc][irad]);
                    for(PS::S64 k = 0; k < NumberOfElement; k++) {
                        fprintf(fp, " %+e", cmps_g[k][iphi][iinc][irad] * ninv);
                        if(k == 3) {
                            fprintf(fp, "\n");
                        }
                    }
                    fprintf(fp, "\n");
                }
            }
        }
        fclose(fp);
#else
        FILE * fp = fopen(ofile, "w");
        PS::S64 iinc = NumberOfInclination/2 - 1;
        PS::F64 inc  = getPositionInclination(iinc);
        for(PS::S64 iphi = 0; iphi < NumberOfAzimuth; iphi++) {
            for(PS::S64 irad = 0; irad < NumberOfRadius; irad++) {
                //PS::F64 rad = getPositionRadius(irad);
                PS::F64 rad = getPositionRadius(irad) / SPHDetailElement::VelocityOverRadius;
                PS::F64 phi = getPositionAzimuth(iphi);
                PS::F64vec pos(rad * sin(inc) * cos(phi),
                               rad * sin(inc) * sin(phi),
                               rad * cos(inc));
                fprintf(fp, "%+e %+e %+e %+e %+e", pos[0], pos[1], pos[2],
                        dens_g[iphi][iinc][irad],
                        temp_g[iphi][iinc][irad]);
                PS::F32 norm = 0.;
                for(PS::S64 k = 0; k < NumberOfElement; k++) {
                    norm += cmps_g[k][iphi][iinc][irad];
                }
                PS::F32 ninv = (norm != 0. && dens_g[iphi][iinc][irad] != 0.) ? (1. / norm) : 0.;
                for(PS::S64 k = 0; k < NumberOfElement; k++) {
                    fprintf(fp, " %+.2e", cmps_g[k][iphi][iinc][irad] * ninv);
                }
                fprintf(fp, "\n");
            }
        }
        fclose(fp);
#endif
    }

}

int main(int argc, char ** argv) {
    PS::Initialize(argc, argv);

    PS::ParticleSystem<SPHAnalysis> sph;
    sph.initialize();
    sph.createParticle(0);
    sph.setNumberOfParticleLocal(0);

    PS::ParticleSystem<SPHDetailElement> sphelem;
    sphelem.initialize();
    sphelem.createParticle(0);
    sphelem.setNumberOfParticleLocal(0);

    char idir[1024], otype[1024], etype[1024];
    PS::S64 itime;
    FILE * fp = fopen(argv[1], "r");
    fscanf(fp, "%s", idir);
    fscanf(fp, "%s", otype);
    fscanf(fp, "%s", etype);
    fscanf(fp, "%lld", &itime);
    fscanf(fp, "%lld", &SPHDetailElement::CompanionCO);
    fscanf(fp, "%lf", &SPHDetailElement::VelocityOverRadius);
    fclose(fp);
    assert(SPHDetailElement::VelocityOverRadius > 0.);

    {        
        char tfile[1024];
        FILE *fp = NULL;
        PS::S64 tdir = 0;
        for(PS::S64 iidir = 0; iidir < 100; iidir++) {
            sprintf(tfile, "%s/t%02d/sph_t%04lld_p%06d_i%06d.dat", idir, iidir, itime,
                    PS::Comm::getNumberOfProc(), 0);
            fp = fopen(tfile, "r");
            if(fp != NULL) {
                tdir = iidir;
                break;
            }
        }
        if(fp == NULL) {
            if(PS::Comm::getRank() == 0) {
                fprintf(stderr, "Not found %s\n", tfile);
            }
            PS::Finalize();
            exit(0);
        }
        fclose(fp);

        char sfile[1024];
        sprintf(sfile, "%s/t%02d/sph_t%04d", idir, tdir, itime);
        sph.readParticleAscii(sfile, "%s_p%06d_i%06d.dat");

        char efile[1024];
        sphelem.readParticleAscii(etype, "%s_p%06d_i%06d.data");

        char ofile[1024];
        sprintf(ofile, "%s.dat", otype);
        assignToMesh(ofile, sph, sphelem);

    }

    PS::Finalize();

    return 0;
}
