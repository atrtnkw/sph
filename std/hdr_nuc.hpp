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
    const PS::S64 NumberOfNucleon = 13;
    const PS::F64 ainv[NumberOfNucleon] = {1/4., 1/12., 1/16., 1/20., 1/24.,
                                           1/28., 1/32., 1/36., 1/40.,
                                           1/44., 1/48., 1/52., 1/56.};
    const PS::F64 zaratio = 0.5;
    const PS::F64 minv[NumberOfNucleon] = {4./4.002603, 12./12.000000, 16./15.994915,
                                           20./19.992440, 24./23.985042, 28./27.976927,
                                           32./31.972071, 36./35.967546, 40./39.962591,
                                           44./43.9596901, 48./47.954032, 52./51.948114,
                                           56./55.942132};    
    PS::F64 enuc    = 0.;
}

namespace NR = NuclearReaction;

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

    static PS::F64 countNumberOfNucleon(PS::F64 * cmps) {
        using namespace NuclearReaction;
        PS::F64 nn = 0;
        for(PS::S32 k = 0; k < NumberOfNucleon; k++) {
            nn += minv[k] * cmps[k];
        }
        return nn;
    }

    static PS::F64 callNuclearReaction(PS::F64 dt,
                                       PS::F64 tt,
                                       PS::F64 dd,
                                       PS::F64 * cmps) {
#if 0
        PS::F64 de = 0.;
        solve_aprox13_(&dt, &tt, &dd, cmps, &de);
        return de;
#else
        PS::F64 n0 = getInstance().countNumberOfNucleon(cmps);
        PS::F64 de = 0.;
        solve_aprox13_(&dt, &tt, &dd, cmps, &de);
        PS::F64 n1 = getInstance().countNumberOfNucleon(cmps);
        de = (1. - n0 / n1) * CodeUnit::SpeedOfLight * CodeUnit::SpeedOfLight;
        //printf("%+e %+e %+e %+e %+e\n", de, n0/n1, n0, n1, CodeUnit::SpeedOfLight);
        return de;
#endif
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

        //PS::F64 de = dt * getInstance().callNuclearReaction(dt, tt, dd, composition);
        PS::F64 de = getInstance().callNuclearReaction(dt, tt, dd, composition);
        PS::F64 denergy = de * UnitOfEnergyInv;

        return denergy;
    }
};

typedef CalcNuclearReactionHydrostatic CalcNRH;
