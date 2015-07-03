#include "eos/WDEOS.hpp"

#include <sys/resource.h>
#include <sys/time.h>
#include <cassert>

double getTime(void){
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 0.001 * 0.001);
}

extern "C" {
    void helmeos2_(double * tt, double * dd, double * abar, double * zbar,
                   double * pp, double * u, double * du, double * cs, bool * eosfail);
}

namespace CodeUnit {
    double SolarRadius        = 6.9599e10; // cm
    double SolarMass          = 1.9891e33; // g
    double UnitOfLength       = SolarRadius * 1.0e-3d;
    double UnitOfMass         = SolarMass   * 1.0e-6d;
    double UnitOfTime         = 1.0d;
    double UnitOfVelocity     = UnitOfLength / UnitOfTime;
    double UnitOfDensity      = UnitOfMass / (UnitOfLength * UnitOfLength * UnitOfLength);
    double UnitOfEnergy       = UnitOfVelocity * UnitOfVelocity;
    double UnitOfPressure     = UnitOfMass * UnitOfLength
        / (UnitOfTime * UnitOfTime * UnitOfLength * UnitOfLength);
    double UnitOfAcceleration = UnitOfVelocity / UnitOfTime;
}

int main(int argc, char **argv)
{
    const int nden = 128;
    double denmin  = 1e-3;
    double denmax  = 1e10;
    double dden    = pow(denmax / denmin, 1. / nden);
    const int ntmp = 6;
    double tmp     = 1e5;

//    double abar = 14.0;
//    double zbar =  7.0;

    {
        double tt, dd, abar, zbar, pp, u, du, cs;
        bool   eosfail;
        dd   = 1e7;
        abar = 12.0;
        zbar = 6.0;
        tt = 1e10;
        helmeos2_(&tt, &dd, &abar, &zbar, &pp, &u, &du, &cs, &eosfail);
        printf("%+e %+e\n", pp, u);
        tt = 1e9;
        helmeos2_(&tt, &dd, &abar, &zbar, &pp, &u, &du, &cs, &eosfail);
        printf("%+e %+e\n", pp, u);
    }

#if 0
    double t1 = getTime();
    {
        double tttol = 1e-5;
        for(double energy  = 1e14; energy < 5e19; energy *= 1e1) {
            char filename[64];
            sprintf(filename, "heos_u%.0e.log", energy);
            FILE *fp = fopen(filename, "w");
            double density = denmin;
            for(int i = 0; i < nden; i++, density *= dden) {
                double ttmax = 1e11;
                double ttmin = 1e5;
                int    cnt = 2;
                double umax, umin;
                double u, pp, du, cs;
                bool   eosfail;            
                helmeos2_(&ttmax, &density, &abar, &zbar, &pp, &umax, &du, &cs, &eosfail);
                helmeos2_(&ttmin, &density, &abar, &zbar, &pp, &umin, &du, &cs, &eosfail);
                
                if(energy > umax) {
                    fprintf(fp, "Too large energy!\n");
                    exit(0);
                } else if (energy < umin) {
                    fprintf(fp, "%+e %+e %+e %+e %+e %4d\n", density, energy, pp, ttmin,
                            (umin - energy) / energy, cnt);
                    continue;
                } 
                
                double ttmid;
                do {
                    ttmid = sqrt(ttmax * ttmin);
                    helmeos2_(&ttmid, &density, &abar, &zbar, &pp, &u, &du, &cs, &eosfail);
                    cnt++;
                    if(u < energy) {
                        ttmin = ttmid;
                    } else {
                        ttmax = ttmid;
                    }
                    assert(cnt < 100);
                }while(ttmax / ttmin > 1. + tttol);
                fprintf(fp, "%+e %+e %+e %+e %+e %4d\n", density, energy, pp, ttmin,
                        (u - energy) / energy, cnt);                    
            }
            fclose(fp);
        }
    }
    printf("%+e\n", getTime() - t1);
#endif

#if 0
    FILE *fp = fopen("heos.log", "w");
    for(int i = 0; i < ntmp; i++, tmp *= 1e1) {
        double den = denmin;
        fprintf(fp, "# temperature: %e\n", tmp);
        for(int j = 0; j < nden; j++, den *= dden) {
            double pp, u, du, cs;
            bool eosfail;
            helmeos2_(&tmp, &den, &abar, &zbar,
                      &pp, &u, &du, &cs, &eosfail);
            fprintf(fp, "%+e %+e %+e\n", den, u, pp);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
#endif

    return 0;
}
