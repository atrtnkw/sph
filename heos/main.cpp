#include "eos/WDEOS.hpp"

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
    
    int    nden    = 128;
    double denmin  = 1e-3;
    double denmax  = 1e10;
    double dden    = pow(denmax / denmin, 1. / nden);
    double enemin  = 1e12;
    double enemax  = 1e30;
    double enetol  = 1e-3;

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

    return 0;
}
