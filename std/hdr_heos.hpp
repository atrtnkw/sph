#pragma once

#include "heos/WDEOS.hpp"

#if 0
extern "C" {
    void helmeos2_(double * tt, double * dd, double * abar, double * zbar,
                   double * pp, double * u, double * du, double * cs, bool * eosfail);
}
#endif

extern "C" {
    void init_flash_helmholtz_(double * coulombin);
    void flash_helmholtz_(double * din,
                          double * ein,
                          double * tin,
                          double * xin,
                          double * pout,
                          double * cout,
                          double * tout,
                          double * sout);
    void flash_helmholtz_e_(double * din,
                            double * tin,
                            double * xin,
                            double * eout);
}

namespace CodeUnit {
    PS::F64 SolarRadius        = 6.9599e10;         // cm
    PS::F64 SolarMass          = 1.9891e33;         // g
    PS::F64 GravityConstant    = 6.6738480e-8;      // [cm^3 g^-1 g^-2]
    PS::F64 SpeedOfLight       = 2.99792458e10;     // cm s^-1
    PS::F64 MassOfNeutron      = 1.67492721184e-24; // g
    PS::F64 MassOfProton       = 1.67262163783e-24; // g
    PS::F64 ElectronVoltToErg  = 1.60217648740e-12;
    PS::F64 MegaElectronVoltToErg  = ElectronVoltToErg * 1e6;
    PS::F64 MegaElectronVoltToGram = MegaElectronVoltToErg / (SpeedOfLight * SpeedOfLight);
    PS::F64 AvogadroConstant   = 6.0221417930e23;
    PS::F64 StefanBoltzmannConstant = 5.670367e-5;
    PS::F64 RadiationConstant       = 4. * StefanBoltzmannConstant / SpeedOfLight;
    PS::F64 RadiationConstantInv    = 1. / RadiationConstant;

    PS::F64 UnitOfLength       = SolarRadius * 1.0e-3;
    PS::F64 UnitOfMass         = SolarMass   * 1.0e-6;
    PS::F64 UnitOfTime         = 1.0;

    PS::F64 UnitOfVelocity     = UnitOfLength / UnitOfTime;
#ifdef FOR_TUBE_TEST
    PS::F64 UnitOfDensity      = UnitOfMass /  UnitOfLength;
    PS::F64 UnitOfEnergy       = UnitOfVelocity * UnitOfVelocity;
    PS::F64 UnitOfPressure     = UnitOfMass * UnitOfLength / (UnitOfTime * UnitOfTime);
#else
    PS::F64 UnitOfDensity      = UnitOfMass / (UnitOfLength * UnitOfLength * UnitOfLength);
    PS::F64 UnitOfEnergy       = UnitOfVelocity * UnitOfVelocity;
    PS::F64 UnitOfPressure     = UnitOfMass * UnitOfLength
        / (UnitOfTime * UnitOfTime * UnitOfLength * UnitOfLength);
#endif
    PS::F64 UnitOfAcceleration = UnitOfVelocity / UnitOfTime;

    PS::F64 UnitOfLengthInv    = 1. / UnitOfLength;
    PS::F64 UnitOfMassInv      = 1. / UnitOfMass;
    PS::F64 UnitOfTimeInv      = 1. / UnitOfTime;
    PS::F64 UnitOfVelocityInv  = 1. / UnitOfVelocity;
    PS::F64 UnitOfDensityInv   = 1. / UnitOfDensity;
    PS::F64 UnitOfEnergyInv    = 1. / UnitOfEnergy;
    PS::F64 UnitOfPressureInv  = 1. / UnitOfPressure;

    PS::F64 MaximumOfTemperature       = 1e11;
    PS::F64 MinimumOfTemperature       = 1e5;
    PS::F64 MaximumOfTemperatureNSE    = 1e10;
    PS::F64 MinimumOfTemperatureNSE    = 5e9;
    PS::F64 BoundaryTemperature        = 1e8;
    PS::F64 TolaranceOfLowTemperature  = 1e-4;
    PS::F64 TolaranceOfHighTemperature = 1e-4;

    // different from OTOO code in time unit !!
    PS::F64 GravityConstantInThisUnit  = GravityConstant /
        ((UnitOfLength * UnitOfLength * UnitOfLength) * UnitOfMassInv
         * (UnitOfTimeInv * UnitOfTimeInv));
    PS::F64 SpeedOfLightInThisUnit = SpeedOfLight * UnitOfVelocityInv;
    PS::F64 MinimumOfDensityExplicitInThisUnit = 3e7 * CodeUnit::UnitOfDensityInv;
    PS::F64 MinimumOfDensityNSEInThisUnit = 2e8 * UnitOfDensityInv;

    // for imbh
    PS::F64 BlackHoleMass;
    PS::F64 BlackHoleMassInThisUnit;

    // for eos
#ifdef COULOMB_CORRECTION
    PS::F64 FractionOfCoulombCorrection = 1.;
#else
    PS::F64 FractionOfCoulombCorrection = 0.;
#endif
}

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

namespace NR = NuclearReaction;

class CalcEquationOfState {
private:
    OTOO::EOS * eos_;
    CalcEquationOfState() {
        using namespace OTOO;
        using namespace CodeUnit;
        if(RP::FlagDamping != 1) {
//            eos_ = new WDEOS_D_E(UnitOfDensity, UnitOfEnergy, UnitOfVelocity,
//                                 UnitOfPressure, 1.0e6);
//            init_flash_helmholtz_();
        } else {
            eos_ = new WDEOSforDumping(UnitOfDensity, UnitOfEnergy, UnitOfVelocity,
                                       UnitOfPressure, 1.0e6);
        }
    };
    ~CalcEquationOfState() {};
    CalcEquationOfState(const CalcEquationOfState & c);    
    CalcEquationOfState & operator=(const CalcEquationOfState & c);
    static CalcEquationOfState & getInstance() {
        static CalcEquationOfState inst;
        return inst;
    }
public:
    static PS::F64 getPressure(PS::F64 density,
                               PS::F64 energy) {
        return getInstance().eos_->GetP(density, energy);
    }
    static PS::F64 getTemperature(PS::F64 density,
                                  PS::F64 energy) {
        return getInstance().eos_->GetT(density, energy);
    }
    static PS::F64 getSoundVelocity(PS::F64 density,
                                    PS::F64 energy) {
        return getInstance().eos_->GetS(density, energy);
    }
    static PS::F64 getEnergy(PS::F64 density,
                             PS::F64 energy) {
        return getInstance().eos_->GetE(density, energy);
    }
    static PS::F64 getEnergyMin(PS::F64 density) {
        return getInstance().eos_->GetEmin2(density);
    }

#if 0
    static PS::F64 getEnergyMin(PS::F64 density,
                                PS::F64 abar,
                                PS::F64 zbar) {
        PS::F64 ttmin = CodeUnit::MinimumOfTemperature;
        PS::F64 dd = density * CodeUnit::UnitOfDensity;
        PS::F64 uu, pp, du, cs;
        bool    eosfail;

        helmeos2_(&ttmin, &dd, &abar, &zbar, &pp, &uu, &du, &cs, &eosfail);        
        PS::F64 eg = uu * CodeUnit::UnitOfEnergyInv;
        
        return eg;
    }
    static PS::F64 getEnergyMax(PS::F64 density,
                                PS::F64 abar,
                                PS::F64 zbar) {
        PS::F64 ttmax = CodeUnit::MaximumOfTemperature;
        PS::F64 dd = density * CodeUnit::UnitOfDensity;
        PS::F64 uu, pp, du, cs;
        bool    eosfail;

        helmeos2_(&ttmax, &dd, &abar, &zbar, &pp, &uu, &du, &cs, &eosfail);
        PS::F64 eg = uu * CodeUnit::UnitOfEnergyInv;
        return eg;
    }
    static void getThermodynamicQuantity(PS::F64 density,
                                         PS::F64 energy,
                                         PS::F64 abar,
                                         PS::F64 zbar,
                                         PS::F64 & pressure,
                                         PS::F64 & soundvelocity,
                                         PS::F64 & temperature,
                                         PS::S64 & counteos) {
        using namespace CodeUnit;

        PS::F64 ttmax = (temperature != 0. && 10. * temperature < MaximumOfTemperature) ?
            10. * temperature : MaximumOfTemperature;
        PS::F64 ttmin = MinimumOfTemperature;
        PS::S32 cnt   = 2;
        PS::F64 umax, umin;
        PS::F64 uu, pp, du, cs;
        bool    eosfail;

        PS::F64 dd = density * UnitOfDensity;
        PS::F64 eg = energy  * UnitOfEnergy;

        helmeos2_(&ttmax, &dd, &abar, &zbar, &pp, &umax, &du, &cs, &eosfail);
        helmeos2_(&ttmin, &dd, &abar, &zbar, &pp, &umin, &du, &cs, &eosfail);

        if(eg > umax) {
            ttmax = MaximumOfTemperature;
            helmeos2_(&ttmax, &dd, &abar, &zbar, &pp, &umax, &du, &cs, &eosfail);
        }

        if(eg > umax) {
            fprintf(stderr, "Too large energy!\n");
            fprintf(stderr, "Density: %+e Energy: %+e Abar: %+e\n", dd, eg, abar);
            PS::Abort();
        } else if (eg > umin) {
            PS::F64 ttmid;
            PS::F64 eginv = 1. / eg;
            PS::F64 tttol;
            PS::F64 tterr;
            do {
                ttmid = sqrt(ttmax * ttmin);
                helmeos2_(&ttmid, &dd, &abar, &zbar, &pp, &uu, &du, &cs, &eosfail);
                cnt++;
                if(uu < eg) {
                    ttmin = ttmid;
                } else {
                    ttmax = ttmid;
                }
                assert(cnt < 100);
                tterr = ttmax / ttmin - 1.;
                tttol = (ttmax < BoundaryTemperature) ?
                    TolaranceOfLowTemperature :
                    TolaranceOfHighTemperature;
            }while(tterr > tttol);
            ttmin = ttmid;
        }
        pressure      = pp * UnitOfPressureInv;
        soundvelocity = cs * UnitOfVelocityInv;
        temperature   = ttmin;
        counteos      = cnt;
    }
#else
    static PS::F64 getEnergyMin(PS::F64 density,
                                NR::Nucleon & composition) {
//        getInstance();
        PS::F64 din  = density * CodeUnit::UnitOfDensity;
        PS::F64 tin  = CodeUnit::MinimumOfTemperature;
        PS::F64 eout;

        flash_helmholtz_e_(&din, &tin, composition.getPointer(), &eout);

        PS::F64 eg = eout * CodeUnit::UnitOfEnergyInv;
        
        return eg;
    }

    static PS::F64 getEnergyMax(PS::F64 density,
                                NR::Nucleon & composition) {
//        getInstance();
        PS::F64 din  = density * CodeUnit::UnitOfDensity;
        PS::F64 tin  = CodeUnit::MaximumOfTemperature;
        PS::F64 eout;

        flash_helmholtz_e_(&din, &tin, composition.getPointer(), &eout);

        PS::F64 eg = eout * CodeUnit::UnitOfEnergyInv;
        
        return eg;
    }

    static PS::F64 getEnergyGivenTemperature(PS::F64 density,
                                             PS::F64 temperature,
                                             NR::Nucleon & composition) {
        PS::F64 din  = density * CodeUnit::UnitOfDensity;
        PS::F64 tin  = temperature;
        PS::F64 eout;

        flash_helmholtz_e_(&din, &tin, composition.getPointer(), &eout);

        PS::F64 eg = eout * CodeUnit::UnitOfEnergyInv;
        
        return eg;
    }

    static void getThermodynamicQuantity(PS::F64 density,
                                         PS::F64 energy,
                                         NR::Nucleon & composition,
                                         PS::F64 & pressure,
                                         PS::F64 & soundvelocity,
                                         PS::F64 & temperature,
                                         PS::F64 & entropy) {
//        getInstance();
        PS::F64 din = density * CodeUnit::UnitOfDensity;
        PS::F64 ein = energy * CodeUnit::UnitOfEnergy;
        PS::F64 tin = (temperature != 0.) ? temperature : 1e9;
        PS::F64 pout, cout, tout, sout;
        flash_helmholtz_(&din, &ein, &tin, composition.getPointer(),
                         &pout, &cout, &tout, &sout);
        pressure      = pout * CodeUnit::UnitOfPressureInv;
        soundvelocity = cout * CodeUnit::UnitOfVelocityInv;
        temperature   = tout;
        entropy       = sout * CodeUnit::UnitOfEnergyInv;
    }

    static void getThermodynamicQuantityRadiationDominant(PS::F64 density,
                                                          PS::F64 energy,
                                                          PS::F64 & pressure,
                                                          PS::F64 & soundvelocity,
                                                          PS::F64 & temperature,
                                                          PS::F64 & entropy) {
        PS::F64 din = density * CodeUnit::UnitOfDensity;
        PS::F64 ein = energy * CodeUnit::UnitOfEnergy;
        PS::F64 tin = (temperature != 0.) ? temperature : 1e9;
        PS::F64 pout, cout, tout, sout;

        pout = 0.3333333333333333 * ein * din;
        tout = pow(3. * pout * CodeUnit::RadiationConstantInv, 0.25);

        pressure      = pout * CodeUnit::UnitOfPressureInv;
        soundvelocity = 1. * CodeUnit::UnitOfVelocityInv;
        temperature   = tout;
        entropy       = 1. * CodeUnit::UnitOfEnergyInv;
    }
#endif
};
