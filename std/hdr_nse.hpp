#pragma once

extern "C" {
    void setup_nse_(void);
    void solve_nse_(double *tin,
                    double *din,
                    double *xout);
}


class CalcNuclearStatisticalEquilibrium {
private:
#if 0
    static const PS::S32 NumberOfNucleon = 47;
#endif
    
    CalcNuclearStatisticalEquilibrium() {
        setup_nse_();
        if(NR::First) {
            using namespace NuclearReaction;
            using namespace CodeUnit;
            for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
                PS::F64 mion = nion[k] * MassOfNeutron + zion[k] * MassOfProton
                    - bion[k] * MegaElectronVoltToGram;
                winv[k]    = 1. / (mion * AvogadroConstant);
                bionerg[k] = AvogadroConstant * (bion[k] * MegaElectronVoltToErg);
            }        
            NR::First = false;
        }
    };
    ~CalcNuclearStatisticalEquilibrium() {};
    CalcNuclearStatisticalEquilibrium(const CalcNuclearStatisticalEquilibrium &c);
    CalcNuclearStatisticalEquilibrium & operator=(const CalcNuclearStatisticalEquilibrium & c);
    static CalcNuclearStatisticalEquilibrium & getInstance() {
        static CalcNuclearStatisticalEquilibrium inst;
        return inst;
    }

#if 0    
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
#endif
    
    static PS::F64 calcBindingEnergy(PS::F64 * cmps) {
        using namespace NuclearReaction;
        PS::F64 e = 0.;
        for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
            PS::F64 ymass = cmps[k] * winv[k];
            e += bionerg[k] * ymass;
        }
        return e;
    }

    static PS::F64 callSolveNse(PS::F64 tt,
                                PS::F64 dd,
                                PS::F64 * cmps) {
#if 0
        PS::F64 xout[getInstance().NumberOfNucleon];
        solve_nse_(&tt, &dd, xout);
        getInstance().convertComposition(xout, composition);
#else
        PS::F64 e0 = getInstance().calcBindingEnergy(cmps);
        solve_nse_(&tt, &dd, cmps);
        PS::F64 e1 = getInstance().calcBindingEnergy(cmps);
        PS::F64 de = e1 - e0;
        return de;
#endif
    }

public:

    static PS::F64 getGeneratedEnergy(PS::F64 density,
                                      PS::F64 temperature,
                                      PS::F64 * composition) {
        using namespace CodeUnit;

        PS::F64 dd = density * UnitOfDensity;;
        PS::F64 tt = temperature;
        PS::F64 de = getInstance().callSolveNse(tt, dd, composition);
        PS::F64 denergy = de * UnitOfEnergyInv;
        return denergy;
    }

#if 0
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
#endif


};

typedef CalcNuclearStatisticalEquilibrium CalcNSE;
