#pragma once

#include "heos/WDEOS.hpp"

extern "C" {
    void helmeos2_(double * tt, double * dd, double * abar, double * zbar,
                   double * pp, double * u, double * du, double * cs, bool * eosfail);
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

    PS::F64 UnitOfLength       = SolarRadius * 1.0e-3d;
    PS::F64 UnitOfMass         = SolarMass   * 1.0e-6d;
    PS::F64 UnitOfTime         = 1.0d;

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

    PS::F64 UnitOfLengthInv    = 1.d / UnitOfLength;
    PS::F64 UnitOfMassInv      = 1.d / UnitOfMass;
    PS::F64 UnitOfTimeInv      = 1.d / UnitOfTime;
    PS::F64 UnitOfVelocityInv  = 1.d / UnitOfVelocity;
    PS::F64 UnitOfDensityInv   = 1.d / UnitOfDensity;
    PS::F64 UnitOfEnergyInv    = 1.d / UnitOfEnergy;
    PS::F64 UnitOfPressureInv  = 1.d / UnitOfPressure;

    //PS::F64 MaximumOfTemperature       = 1e10;
    PS::F64 MaximumOfTemperature       = 1e11;
    PS::F64 MinimumOfTemperature       = 1e5;
    PS::F64 BoundaryTemperature        = 1e8;
    PS::F64 TolaranceOfLowTemperature  = 1e-4;
    PS::F64 TolaranceOfHighTemperature = 1e-4;

    // different from OTOO code in time unit !!
    PS::F64 GravityConstantInThisUnit  = GravityConstant /
        ((UnitOfLength * UnitOfLength * UnitOfLength) * UnitOfMassInv
         * (UnitOfTimeInv * UnitOfTimeInv));
    PS::F64 SpeedOfLightInThisUnit = SpeedOfLight * UnitOfVelocityInv;
}

class CalcEquationOfState {
private:
    OTOO::EOS * eos_;
    CalcEquationOfState() {
        using namespace OTOO;
        using namespace CodeUnit;
        if(RP::FlagDamping != 1) {
            eos_ = new WDEOS_D_E(UnitOfDensity, UnitOfEnergy, UnitOfVelocity,
                                 UnitOfPressure, 1.0e6);
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

        PS::F64 ttmax = (temperature != 0.d && 10.d * temperature < MaximumOfTemperature) ?
            10.d * temperature : MaximumOfTemperature;
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
            PS::F64 eginv = 1.d / eg;
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
                tterr = ttmax / ttmin - 1.d;
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

};
