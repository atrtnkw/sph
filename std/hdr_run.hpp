#pragma once

namespace RunParameter {
    PS::F64    Time;
    PS::F64    TimeEnd;
    PS::F64    Timestep;
    PS::F64    MaximumTimestep;
    PS::F64    MinimumTimestep;
    PS::F64    TimeAscii;
    PS::F64    TimestepAscii;
    PS::F64    TimeHexa;
    PS::F64    TimestepHexa;
    PS::S64    NumberOfStep;
    PS::S64    NumberOfAscii;
    PS::S64    NumberOfHexa;
    PS::F64ort ComputationalBox;
    PS::S64    NumberOfDimension;
    KernelType KernelType;

    PS::F64    CoefficientOfTimestep;
    PS::F64    GravitationalSoftening;
    PS::F64    AlphaMaximum, AlphaMinimum, AlphaInit;
    PS::F64    AlphuMaximum, AlphuMinimum, AlphuInit;
    PS::F64    KernelSupportRadiusMaximum; // ksrmax    
    PS::F64    EpsilonOfInternalEnergy;    // epsu
    PS::F64    ReductionTime;
    PS::F64    ReductionTimeInv;
    
    PS::S64    NumberOfParticle;
    PS::F64    AdiabaticIndex;     // gamma
    PS::F64vec RotationalVelocity; // omega
    PS::F64    NuclearEnergyTotal; // enuc

    PS::S64    FlagGravity; // 0: off, 1: on
    PS::S64    FlagDamping; // 0: normal, 1: single, 2: binary
    PS::S64    FlagNuclear; // 0: off, 1: on

    FILE * FilePointerForLog;
    FILE * FilePointerForTime;

    PS::S64 readAscii(const char * const filename);
    template <class Tdinfo>
    void readHexa(const char * const filename,
                  Tdinfo & dinfo);
    template <class Tdinfo>
    void writeHexa(const char * const filename,
                   Tdinfo & dinfo);
};

namespace RP = RunParameter;
