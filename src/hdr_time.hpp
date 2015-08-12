#include <sys/stat.h>
#include <sys/time.h>

PS::F64 getWallclockTime() {
    struct timeval tv;
    gettimeofday(& tv, NULL);
    return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6);
}
    
class WallclockTime {
    PS::F64 TimeDecomposeDomain;
    PS::F64 TimeExchangeParticle;
    PS::F64 TimeCalcDensity;
    PS::F64 TimeCalcHydro;
    PS::F64 TimeCalcGravity;
    PS::F64 TimeReferEquationOfState;
    PS::F64 TimeIntegrateOrbit;
    PS::F64 TimeOthers;
    PS::F64 tstart;

    WallclockTime() {}
    ~WallclockTime() {}
    WallclockTime(const WallclockTime & c);
    WallclockTime & operator = (const WallclockTime & c);
    static WallclockTime & getInstance() {
        static WallclockTime inst;
        return inst;
    }

    PS::F64 getWallclockTime() {
        struct timeval tv;
        gettimeofday(& tv, NULL);
        return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6);
    }

public:
    static void start() {
        getInstance().tstart = getInstance().getWallclockTime();
    }
    
    static PS::F64 getTimeTotal() {
        WallclockTime & p = getInstance();
        PS::F64 ttotal;
        ttotal  = p.TimeDecomposeDomain;
        ttotal += p.TimeExchangeParticle;
        ttotal += p.TimeCalcDensity;
        ttotal += p.TimeCalcHydro;
        ttotal += p.TimeCalcGravity;
        ttotal += p.TimeReferEquationOfState;
        ttotal += p.TimeIntegrateOrbit;
        ttotal += p.TimeOthers;
        return ttotal;
    }
    
    static void accumulateDecomposeDomain() {
        WallclockTime & p = getInstance();
        p.TimeDecomposeDomain += p.getWallclockTime() - p.tstart;
    }

    static void accumulateExchangeParticle() {
        WallclockTime & p = getInstance();
        p.TimeExchangeParticle += p.getWallclockTime() - p.tstart;
    }

    static void accumulateCalcDensity() {
        WallclockTime & p = getInstance();
        p.TimeCalcDensity += p.getWallclockTime() - p.tstart;
    }

    static void accumulateCalcHydro() {
        WallclockTime & p = getInstance();
        p.TimeCalcHydro += p.getWallclockTime() - p.tstart;
    }

    static void accumulateCalcGravity() {
        WallclockTime & p = getInstance();
        p.TimeCalcGravity += p.getWallclockTime() - p.tstart;
    }

    static void accumulateReferEquationOfState() {
        WallclockTime & p = getInstance();
        p.TimeReferEquationOfState += p.getWallclockTime() - p.tstart;
    }

    static void accumulateIntegrateOrbit() {
        WallclockTime & p = getInstance();
        p.TimeIntegrateOrbit += p.getWallclockTime() - p.tstart;
    }

    static void accumulateOthers() {
        WallclockTime & p = getInstance();
        p.TimeOthers += p.getWallclockTime() - p.tstart;
    }

    static void reduceInterProcess() {
        WallclockTime & p = getInstance();
        PS::F64 MaxTimeDecomposeDomain       = PS::Comm::getMaxValue(p.TimeDecomposeDomain);
        PS::F64 MaxTimeExchangeParticle      = PS::Comm::getMaxValue(p.TimeExchangeParticle);
        PS::F64 MaxTimeCalcDensity           = PS::Comm::getMaxValue(p.TimeCalcDensity);
        PS::F64 MaxTimeCalcHydro             = PS::Comm::getMaxValue(p.TimeCalcHydro);
        PS::F64 MaxTimeCalcGravity           = PS::Comm::getMaxValue(p.TimeCalcGravity);
        PS::F64 MaxTimeReferEquationOfState  = PS::Comm::getMaxValue(p.TimeReferEquationOfState);
        PS::F64 MaxTimeIntegrateOrbit        = PS::Comm::getMaxValue(p.TimeIntegrateOrbit);
        PS::F64 MaxTimeOthers                = PS::Comm::getMaxValue(p.TimeOthers);
        p.TimeDecomposeDomain      = MaxTimeDecomposeDomain;
        p.TimeExchangeParticle     = MaxTimeExchangeParticle;
        p.TimeCalcDensity          = MaxTimeCalcDensity;
        p.TimeCalcHydro            = MaxTimeCalcHydro;
        p.TimeCalcGravity          = MaxTimeCalcGravity;
        p.TimeReferEquationOfState = MaxTimeReferEquationOfState;
        p.TimeIntegrateOrbit       = MaxTimeIntegrateOrbit;
        p.TimeOthers               = MaxTimeOthers;
    }

    static void dump(PS::F64 time, FILE * fp) {
        WallclockTime & p = getInstance();
        fprintf(fp, "\n");
        fprintf(fp, "Time: %e\n", time);
        fprintf(fp, "TimeDecomposeDomain:          %e\n", p.TimeDecomposeDomain);
        fprintf(fp, "TimeExchangeParticle:         %e\n", p.TimeExchangeParticle);
        fprintf(fp, "TimeCalcDensity:              %e\n", p.TimeCalcDensity);
        fprintf(fp, "TimeCalcHydro:                %e\n", p.TimeCalcHydro);
        fprintf(fp, "TimeCalcGravity:              %e\n", p.TimeCalcGravity);
        fprintf(fp, "TimeCalcReferEquationOfState: %e\n", p.TimeReferEquationOfState);
        fprintf(fp, "TimeIntegrationOrbit:         %e\n", p.TimeIntegrateOrbit);
        fprintf(fp, "TimeOthers:                   %e\n", p.TimeOthers);
        fprintf(fp, "TimeTotal:                    %e\n", p.getTimeTotal());
    }

    static void clear() {
        WallclockTime & p = getInstance();
        p.TimeDecomposeDomain      = 0.d;
        p.TimeExchangeParticle     = 0.d;
        p.TimeCalcDensity          = 0.d;
        p.TimeCalcHydro            = 0.d;
        p.TimeCalcGravity          = 0.d;
        p.TimeReferEquationOfState = 0.d;
        p.TimeIntegrateOrbit       = 0.d;
        p.TimeOthers               = 0.d;
    }
};

typedef WallclockTime WT;
