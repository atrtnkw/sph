#pragma once

#include "heos/WDEOS.hpp"

namespace CodeUnit {
    PS::F64 SolarRadius        = 6.9599e10; // cm
    PS::F64 SolarMass          = 1.9891e33; // g
    PS::F64 GravityConstant    = 6.6738480e-8;  // [cm^3 g^-1 g^-2]

    PS::F64 UnitOfLength       = SolarRadius * 1.0e-3d;
    PS::F64 UnitOfMass         = SolarMass   * 1.0e-6d;
    PS::F64 UnitOfTime         = 1.0d;

    PS::F64 UnitOfVelocity     = UnitOfLength / UnitOfTime;
    PS::F64 UnitOfDensity      = UnitOfMass / (UnitOfLength * UnitOfLength * UnitOfLength);
    PS::F64 UnitOfEnergy       = UnitOfVelocity * UnitOfVelocity;
    PS::F64 UnitOfPressure     = UnitOfMass * UnitOfLength
        / (UnitOfTime * UnitOfTime * UnitOfLength * UnitOfLength);
    PS::F64 UnitOfAcceleration = UnitOfVelocity / UnitOfTime;

    PS::F64 UnitOfLengthInv    = 1.d / UnitOfLength;
    PS::F64 UnitOfMassInv      = 1.d / UnitOfMass;
    PS::F64 UnitOfTimeInv      = 1.d / UnitOfTime;
    PS::F64 UnitOfVelocityInv  = 1.d / UnitOfVelocity;
    PS::F64 UnitOfDensityInv   = 1.d / UnitOfDensity;
    PS::F64 UnitOfEnergyInv    = 1.d / UnitOfEnergy;
    PS::F64 UnitOfPressureInv  = 1.d / UnitOfPressure;

    // different from OTOO code in time unit !!
    PS::F64 grav               = GravityConstant /
        ((UnitOfLength * UnitOfLength * UnitOfLength) * UnitOfMassInv
         * (UnitOfTimeInv * UnitOfTimeInv));
}

class CalcEquationOfState {
private:
    OTOO::EOS * eos_;
    CalcEquationOfState() {
        using namespace OTOO;
        using namespace CodeUnit;
#ifdef WD_DAMPING1
        eos_ = new WDEOSforDumping(UnitOfDensity, UnitOfEnergy, UnitOfVelocity,
                                   UnitOfPressure, 1.0e6);
#else
        eos_ = new WDEOS_D_E(UnitOfDensity, UnitOfEnergy, UnitOfVelocity,
                             UnitOfPressure, 1.0e6);
#endif
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
};
