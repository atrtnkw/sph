#pragma once

extern "C" {
    void setup_nse_(void);
    void solve_nse_(double *tin,
                    double *din,
                    double *xout);
}


class CalcNuclearStatisticalEquilibrium {
private:
    static const PS::S32 NumberOfNucleon = 47;
    
    CalcNuclearStatisticalEquilibrium() {
        setup_nse_();
    };
    ~CalcNuclearStatisticalEquilibrium() {};
    CalcNuclearStatisticalEquilibrium(const CalcNuclearStatisticalEquilibrium &c);
    CalcNuclearStatisticalEquilibrium & operator=(const CalcNuclearStatisticalEquilibrium & c);
    static CalcNuclearStatisticalEquilibrium & getInstance() {
        static CalcNuclearStatisticalEquilibrium inst;
        return inst;
    }
    
    static void convertComposition(PS::F64 * xout,
                                   PS::F64 * composition) {
        composition[0]  = xout[1];
        composition[1]  = xout[2];
        composition[2]  = xout[8];
        composition[3]  = xout[16];
        composition[4]  = xout[18];
        composition[5]  = xout[20];
        composition[6]  = xout[23];
        composition[7]  = xout[25];
        composition[8]  = xout[27];
        composition[9]  = xout[28];
        composition[10] = xout[29];
        composition[11] = xout[36];
        composition[12] = xout[43];
        PS::F64 norm = 0.;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            norm += composition[k];
        }
        PS::F64 norminv = 1. / norm;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            composition[k] *= norminv;
        }
    }
    
    static void callSolveNse(PS::F64 tt,
                             PS::F64 dd,
                             PS::F64 * composition) {
        PS::F64 xout[getInstance().NumberOfNucleon];
        solve_nse_(&tt, &dd, xout);
        getInstance().convertComposition(xout, composition);
    }

public:

    static PS::F64 getGeneratedEnergy(PS::F64 density,
                                      PS::F64 temperature,
                                      PS::F64 * composition) {
        using namespace CodeUnit;

        PS::F64 dd = density * UnitOfDensity;;
        PS::F64 tt = temperature;
        getInstance().callSolveNse(tt, dd, composition);
        return 0.;
    }

    static bool IsNseConditionSatisfied(PS::F64 density,
                                        PS::F64 temperature) {
        using namespace CodeUnit;

        bool x;
        PS::F64 dd = density * UnitOfDensity;;
        PS::F64 tt = temperature;
        if(dd > 5e7 && tt > 6e9) {
            x = true;
        } else {
            x = false;
        }
        return x;
    }


};

typedef CalcNuclearStatisticalEquilibrium CalcNSE;
