#pragma once

namespace CodeUnit {
    PS::F64 SolarRadius        = 6.9599e10; // cm
    PS::F64 SolarMass          = 1.9891e33; // g
    PS::F64 GravityConstant    = 6.6738480e-8;  // [cm^3 g^-1 g^-2]

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

    PS::F64 MaximumOfTemperature       = 1e10;
    PS::F64 MinimumOfTemperature       = 1e5;
    PS::F64 BoundaryTemperature        = 1e8;
    PS::F64 TolaranceOfLowTemperature  = 1e-4;
    PS::F64 TolaranceOfHighTemperature = 1e-4;

    // different from OTOO code in time unit !!
    PS::F64 GravityConstantInThisUnit  = GravityConstant /
        ((UnitOfLength * UnitOfLength * UnitOfLength) * UnitOfMassInv
         * (UnitOfTimeInv * UnitOfTimeInv));
}
