#pragma once

extern "C" {
    void setup_aprox13_(void);
    void solve_aprox13_(double *tstep,
                        double *tin,
                        double *din,
                        double *xin,
                        double *deout);
}

namespace NuclearReaction {
    static const PS::S32 NumberOfNucleon = 13;
}

class CalcNuclearReactionHydrostatic {
private:
    CalcNuclearReactionHydrostatic() {
        setup_aprox13_();
    };
    ~CalcNuclearReactionHydrostatic() {};
    CalcNuclearReactionHydrostatic(const CalcNuclearReactionHydrostatic &c);
    CalcNuclearReactionHydrostatic & operator=(const CalcNuclearReactionHydrostatic & c);
    static CalcNuclearReactionHydrostatic & getInstance() {
        static CalcNuclearReactionHydrostatic inst;
        return inst;
    }
    static PS::F64 callNuclearReaction(PS::F64 dt,
                                       PS::F64 tt,
                                       PS::F64 dd,
                                       PS::F64 * cmps) {
        PS::F64 de = 0.;
        solve_aprox13_(&dt, &tt, &dd, cmps, &de);
        return de;
    }
public:
    static PS::F64 getGeneratedEnergy(PS::F64 dtime,
                                      PS::F64 density,
                                      PS::F64 temperature,
                                      PS::F64 * composition) {
        using namespace CodeUnit;
        using namespace NuclearReaction;

        PS::F64 dt = dtime * UnitOfTimeInv;
        PS::F64 dd = density * UnitOfDensity;;
        PS::F64 tt = temperature;

        PS::F64 de = dt * getInstance().callNuclearReaction(dt, tt, dd, composition);
        PS::F64 denergy = de * UnitOfEnergyInv;

        return denergy;
    }
};

typedef CalcNuclearReactionHydrostatic CalcNRH;
