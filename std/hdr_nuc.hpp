#pragma once

extern "C" {
    void setup_aprox13_(void);
    void solve_aprox13_(double *tstep,
                        double *tin,
                        double *din,
                        double *xin,
                        double *deout);
}

extern "C" {
    void init_flash_burner_();
    void bn_burner_(double * tstep,
                    double * tin,
                    double * din,
                    double * xin,
                    double * xout,
                    double * sdot);
}

/*
namespace NuclearReaction {    
    const PS::S64 NumberOfNucleon = 13;
    const PS::F64 ainv[NumberOfNucleon] = {1/4., 1/12., 1/16., 1/20., 1/24.,
                                           1/28., 1/32., 1/36., 1/40.,
                                           1/44., 1/48., 1/52., 1/56.};
    const PS::F64 zion[NumberOfNucleon] = {2.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0,
                                           18.0, 20.0, 22.0, 24.0, 26.0, 28.0};
    const PS::F64 nion[NumberOfNucleon] = {2.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0,
                                           18.0, 20.0, 22.0, 24.0, 26.0, 28.0};
    const PS::F64 bion[NumberOfNucleon] = {28.29603, 92.16294, 127.62093, 160.64788,
                                           198.25790, 236.53790, 271.78250, 306.72020,
                                           342.05680, 375.47720, 411.46900, 447.70800,
                                           484.00300};
    const PS::F64 zaratio = 0.5;
    const PS::F64 minv[NumberOfNucleon] = {4./4.002603, 12./12.000000, 16./15.994915,
                                           20./19.992440, 24./23.985042, 28./27.976927,
                                           32./31.972071, 36./35.967546, 40./39.962591,
                                           44./43.9596901, 48./47.954032, 52./51.948114,
                                           56./55.942132};    
    PS::F64 enuc = 0.;
    PS::F64 winv[NumberOfNucleon]    = {0., 0., 0., 0., 0., 0., 0.,
                                        0., 0., 0., 0., 0., 0.};
    PS::F64 bionerg[NumberOfNucleon] = {0., 0., 0., 0., 0., 0., 0.,
                                        0., 0., 0., 0., 0., 0.};;
    bool First = true;

    class Nucleon {
        PS::F64 val[NumberOfNucleon];
    public:
        Nucleon () {
            for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
                val[k] = 0.;
            }
        }
        Nucleon (const Nucleon & src) {
            for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
                val[k] = src.val[k];
            }
        }
        PS::F64 & operator [] (const PS::S32 k) {return val[k];}
        const PS::F64 & operator [] (const PS::S32 k) const {return val[k];}
        const Nucleon & operator = (const Nucleon & src) {
            for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
                val[k] = src.val[k];
            }
            return (*this);
        }
        PS::F64 * getPointer() {return val;}
        void print(const PS::S64 id = 0) {
            printf("%10d", id);
            for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
                printf(" %.3e", val[k]);
            }
            printf("\n");
        }
    };
}
*/

namespace NR = NuclearReaction;

class CalcNuclearReactionHydrostatic {
private:
    CalcNuclearReactionHydrostatic() {
        /*
        setup_aprox13_();
        if(NR::First) {
            using namespace NuclearReaction;
            using namespace CodeUnit;
            for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
                PS::F64 mion = nion[k] * MassOfNeutron + zion[k] * MassOfProton
                    - bion[k] * MegaElectronVoltToGram;
                //winv[k]    = 1. / (mion * AvogadroConstant);
                winv[k]    = ainv[k];
                bionerg[k] = AvogadroConstant * (bion[k] * MegaElectronVoltToErg);
            }        
            NR::First = false;
        }
        */
        init_flash_burner_();
    };
    ~CalcNuclearReactionHydrostatic() {};
    CalcNuclearReactionHydrostatic(const CalcNuclearReactionHydrostatic &c);
    CalcNuclearReactionHydrostatic & operator=(const CalcNuclearReactionHydrostatic & c);
    static CalcNuclearReactionHydrostatic & getInstance() {
        static CalcNuclearReactionHydrostatic inst;
        return inst;
    }

    /*
    static PS::F64 calcBindingEnergy(PS::F64 * cmps) {
        using namespace NuclearReaction;
        PS::F64 e = 0.;
        for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
            PS::F64 ymass = cmps[k] * winv[k];
            e += bionerg[k] * ymass;
        }
        return e;
    }

    static PS::F64 callNuclearReaction(PS::F64 dt,
                                       PS::F64 tt,
                                       PS::F64 dd,
                                       PS::F64 * cmps) {
        PS::F64 e0  = getInstance().calcBindingEnergy(cmps);
        PS::F64 dum = 0.;
        solve_aprox13_(&dt, &tt, &dd, cmps, &dum);
        PS::F64 e1  = getInstance().calcBindingEnergy(cmps);
        PS::F64 de = e1 - e0;
        return de;
    }
    */

    static PS::F64 callNuclearReaction(PS::F64 dt,
                                       PS::F64 tt,
                                       PS::F64 dd,
                                       NR::Nucleon & cmps0) {
        NR::Nucleon cmps1;
        PS::F64 sdot;
        bn_burner_(&dt, &tt, &dd, cmps0.getPointer(), cmps1.getPointer(), &sdot);
        PS::F64 de = dt * sdot;
        cmps0 = cmps1;
        return de;
    }

public:
    /*
    static PS::F64 getGeneratedEnergy(PS::F64 dtime,
                                      PS::F64 density,
                                      PS::F64 temperature,
                                      PS::F64 * composition) {
        using namespace CodeUnit;

        PS::F64 dt = dtime * UnitOfTime;
        PS::F64 dd = density * UnitOfDensity;;
        PS::F64 tt = temperature;
        PS::F64 de = getInstance().callNuclearReaction(dt, tt, dd, composition);
        PS::F64 denergy = de * UnitOfEnergyInv;
        return denergy;
    }
    */
    static PS::F64 getGeneratedEnergy(PS::F64 dtime,
                                      PS::F64 density,
                                      PS::F64 temperature,
                                      NR::Nucleon & composition) {
        using namespace CodeUnit;

        PS::F64 dt = dtime * UnitOfTime;
        PS::F64 dd = density * UnitOfDensity;;
        PS::F64 tt = temperature;
        PS::F64 de = getInstance().callNuclearReaction(dt, tt, dd, composition);
        PS::F64 denergy = de * UnitOfEnergyInv;
        return denergy;
    }
};

typedef CalcNuclearReactionHydrostatic CalcNRH;
