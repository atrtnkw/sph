#pragma once

extern "C" {
    void setup_iso7_(void);
    void solve_iso7_(double *tstep,
                     double *tin,
                     double *din,
                     double *xin,
                     double *deout);
}

namespace NuclearReaction {
    const PS::S64 NumberOfNucleon = 7;
    const PS::F64 ainv[NumberOfNucleon] = {1/4., 1/12., 1/16., 1/20., 1/24.,
                                           1/28., 1/56.};
    const PS::F64 zaratio = 0.5;
    PS::F64 enuc    = 0.;
}

namespace NR = NuclearReaction;

class CalcNuclearReactionHydrostatic {
private:
    CalcNuclearReactionHydrostatic() {
        setup_iso7_();
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
        solve_iso7_(&dt, &tt, &dd, cmps, &de);
        return de;
    }
public:
    static PS::F64 getGeneratedEnergy(PS::F64 dtime,
                                      PS::F64 density,
                                      PS::F64 temperature,
                                      PS::F64 * composition) {
        using namespace CodeUnit;
        using namespace NuclearReaction;

        PS::F64 dt = dtime * UnitOfTime;
        PS::F64 dd = density * UnitOfDensity;;
        PS::F64 tt = temperature;

        PS::F64 de = getInstance().callNuclearReaction(dt, tt, dd, composition);
        PS::F64 denergy = de * UnitOfEnergyInv;

        return denergy;
    }
};

typedef CalcNuclearReactionHydrostatic CalcNRH;
