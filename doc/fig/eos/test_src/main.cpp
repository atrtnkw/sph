#include "eos/WDEOS.hpp"

#include <sys/resource.h>
#include <sys/time.h>

double getTime(void){
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 0.001 * 0.001);
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
    using namespace OTOO;
    using namespace CodeUnit;

    EOS *eos;

    eos = new WDEOS_D_E(UnitOfDensity, UnitOfEnergy, UnitOfVelocity,
                        UnitOfPressure, 1.0e6);
    
//    const int nden = 128;
    const int nden = 1048576;
    double denmin  = 1e-3;
    double denmax  = 1e10;
    double dden    = pow(denmax / denmin, 1. / nden);
    double enemin  = 1e12;
    double enemax  = 1e30;
    double enetol  = 1e-3;

    {
        static double d[nden], p[nden], u[nden], c[nden];
        double dd = denmin;
        double tt = 1e7;
        double t1 = getTime();
        for(int i = 0; i < nden; i++, dd *= dden) {
            d[i] = dd;
            eosx_return_(&tt, &dd, &p[i], &u[i], &c[i]);
        }
        fprintf(stderr, "%+e\n", getTime() - t1);
        for(int i = 0; i < nden; i++) {
            if(i % 100 == 0)
                printf("%+e %+e %+e %+e %+e\n", tt, d[i], p[i], u[i], c[i]);
        }        
    }

#if 0
    FILE *fp = fopen("heos.log", "w");
    int    ntemp = 6;
    double temp  = 1e5;
    for(int itemp = 0; itemp < ntemp; itemp++, temp *= 1e1) {
        double den     = denmin;
        for(int i = 0; i < nden; i++, den *= dden) {
            double em, tm;
            double el = enemin;
            double eu = enemax / (den / denmin);

            eu = (eu > 1e19) ? eu : 1e19;
            do {
                em = sqrt(el * eu);
                tm = eos->GetT(den / UnitOfDensity, em / UnitOfEnergy);
                if(temp < tm) {
                    eu = em;
                } else {
                    el = em;
                }
            }while(eu / el > 1. + enetol);

            fprintf(fp, "%+e %+e %+e %+e\n", den, em,
                   eos->GetP(den / UnitOfDensity, em / UnitOfEnergy) * UnitOfPressure,
                   temp);
            
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
#endif

    return 0;
}
