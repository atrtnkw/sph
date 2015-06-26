#include "eos/WDEOS.hpp"

#include <sys/resource.h>
#include <sys/time.h>

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

    double abar = 56.0;
    double zbar = 28.0;

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

    return 0;
}
