#pragma once

namespace RunParameter {
    PS::F64    Time;
    PS::F64    TimeEnd;
    PS::F64    Timestep;
    PS::F64    TimestepAscii;
    PS::F64    TimestepHexa;
    PS::S64    NumberOfStep;
    PS::F64ort ComputationalBox;
    PS::S64    NumberOfDimension;
    KernelType KernelType;

    PS::F64    CoefficientOfTimestep;
    PS::F64    GravitationalSoftening;
    PS::F64    AlphaMaximum, AlphaMinimum;
    PS::F64    KernelSupportRadiusMaximum; // ksrmax    
    
    PS::S64    NumberOfParticle;
    PS::F64    AdiabaticIndex;     // gamma
    PS::F64vec RotationalVelocity; // omega
    PS::F64    NuclearEnergyTotal; // enuc

    PS::S64 readAscii(const char * const filename);
    template <class Tdinfo>
    void readHexa(const char * const filename,
                  Tdinfo & dinfo);
    template <class Tdinfo>
    void writeHexa(const char * const filename,
                   Tdinfo & dinfo);
};

namespace RP = RunParameter;
