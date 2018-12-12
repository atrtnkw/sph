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
    PS::F64 TimeCalcNuclearReaction;
    PS::F64 TimeOthers;
    PS::F64 tstart;

    PS::F64 TimeDecomposeDomainAverage;
    PS::F64 TimeExchangeParticleAverage;
    PS::F64 TimeCalcDensityAverage;
    PS::F64 TimeCalcHydroAverage;
    PS::F64 TimeCalcGravityAverage;
    PS::F64 TimeReferEquationOfStateAverage;
    PS::F64 TimeIntegrateOrbitAverage;
    PS::F64 TimeCalcNuclearReactionAverage;
    PS::F64 TimeOthersAverage;

    PS::F64 TimeDecomposeDomainMax;
    PS::F64 TimeExchangeParticleMax;
    PS::F64 TimeCalcDensityMax;
    PS::F64 TimeCalcHydroMax;
    PS::F64 TimeCalcGravityMax;
    PS::F64 TimeReferEquationOfStateMax;
    PS::F64 TimeIntegrateOrbitMax;
    PS::F64 TimeCalcNuclearReactionMax;
    PS::F64 TimeOthersMax;

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
#if 0
        ttotal  = p.TimeDecomposeDomain;
        ttotal += p.TimeExchangeParticle;
        ttotal += p.TimeCalcDensity;
        ttotal += p.TimeCalcHydro;
        ttotal += p.TimeCalcGravity;
        ttotal += p.TimeReferEquationOfState;
        ttotal += p.TimeIntegrateOrbit;
        ttotal += p.TimeCalcNuclearReaction;
        ttotal += p.TimeOthers;
#else
        ttotal  = p.TimeDecomposeDomainMax;
        ttotal += p.TimeExchangeParticleMax;
        ttotal += p.TimeCalcDensityMax;
        ttotal += p.TimeCalcHydroMax;
        ttotal += p.TimeCalcGravityMax;
        ttotal += p.TimeReferEquationOfStateMax;
        ttotal += p.TimeIntegrateOrbitMax;
        ttotal += p.TimeCalcNuclearReactionMax;
        ttotal += p.TimeOthersMax;
#endif
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

    static void accumulateCalcNuclearReaction() {
        WallclockTime & p = getInstance();
        p.TimeCalcNuclearReaction += p.getWallclockTime() - p.tstart;
    }

    static void accumulateOthers() {
        WallclockTime & p = getInstance();
        p.TimeOthers += p.getWallclockTime() - p.tstart;
    }

    static void reduceInterProcess() {
        WallclockTime & p = getInstance();
        p.TimeDecomposeDomainMax       = PS::Comm::getMaxValue(p.TimeDecomposeDomain);
        p.TimeExchangeParticleMax      = PS::Comm::getMaxValue(p.TimeExchangeParticle);
        p.TimeCalcDensityMax           = PS::Comm::getMaxValue(p.TimeCalcDensity);
        p.TimeCalcHydroMax             = PS::Comm::getMaxValue(p.TimeCalcHydro);
        p.TimeCalcGravityMax           = PS::Comm::getMaxValue(p.TimeCalcGravity);
        p.TimeReferEquationOfStateMax  = PS::Comm::getMaxValue(p.TimeReferEquationOfState);
        p.TimeIntegrateOrbitMax        = PS::Comm::getMaxValue(p.TimeIntegrateOrbit);
        p.TimeCalcNuclearReactionMax   = PS::Comm::getMaxValue(p.TimeCalcNuclearReaction);
        p.TimeOthersMax                = PS::Comm::getMaxValue(p.TimeOthers);

        PS::F64 rankinv = 1.0 / (PS::F64)(PS::Comm::getNumberOfProc());
        p.TimeDecomposeDomainAverage      = PS::Comm::getSum(p.TimeDecomposeDomain) * rankinv;
        p.TimeExchangeParticleAverage     = PS::Comm::getSum(p.TimeExchangeParticle) * rankinv;
        p.TimeCalcDensityAverage          = PS::Comm::getSum(p.TimeCalcDensity) * rankinv;
        p.TimeCalcHydroAverage            = PS::Comm::getSum(p.TimeCalcHydro) * rankinv;
        p.TimeCalcGravityAverage          = PS::Comm::getSum(p.TimeCalcGravity) * rankinv;
        p.TimeReferEquationOfStateAverage = PS::Comm::getSum(p.TimeReferEquationOfState) * rankinv;
        p.TimeIntegrateOrbitAverage       = PS::Comm::getSum(p.TimeIntegrateOrbit) * rankinv;
        p.TimeCalcNuclearReactionAverage  = PS::Comm::getSum(p.TimeCalcNuclearReaction) * rankinv;
        p.TimeOthersAverage               = PS::Comm::getSum(p.TimeOthers) * rankinv;        
    }

    static void dump(PS::F64 time, FILE * fp) {
        static bool first = true;
        WallclockTime & p = getInstance();
        if(first) {
            fprintf(fp, "# Time,");
            fprintf(fp, " TimeTotal");
            fprintf(fp, " TimeDecomposeDomain");
            fprintf(fp, " TimeExchangeParticle");
            fprintf(fp, " TimeCalcDensity");
            fprintf(fp, " TimeCalcHydro");
            fprintf(fp, " TimeCalcGravity");
            fprintf(fp, " TimeCalcReferEquationOfState");
            fprintf(fp, " TimeIntegrationOrbit");
            fprintf(fp, " TimeCalcNuclearReaction");
            fprintf(fp, " TimeOthers");
            fprintf(fp, "\n");
            first = false;
        }
        fprintf(fp, " %e", time);
        fprintf(fp, " %e", p.getTimeTotal());
        fprintf(fp, " %e", p.TimeDecomposeDomainMax);
        fprintf(fp, " %e", p.TimeExchangeParticleMax);
        fprintf(fp, " %e", p.TimeCalcDensityMax);
        fprintf(fp, " %e", p.TimeCalcHydroMax);
        fprintf(fp, " %e", p.TimeCalcGravityMax);
        fprintf(fp, " %e", p.TimeReferEquationOfStateMax);
        fprintf(fp, " %e", p.TimeIntegrateOrbitMax);
        fprintf(fp, " %e", p.TimeCalcNuclearReactionMax);
        fprintf(fp, " %e", p.TimeOthersMax);
        fprintf(fp, " %e", p.TimeDecomposeDomainAverage);
        fprintf(fp, " %e", p.TimeExchangeParticleAverage);
        fprintf(fp, " %e", p.TimeCalcDensityAverage);
        fprintf(fp, " %e", p.TimeCalcHydroAverage);
        fprintf(fp, " %e", p.TimeCalcGravityAverage);
        fprintf(fp, " %e", p.TimeReferEquationOfStateAverage);
        fprintf(fp, " %e", p.TimeIntegrateOrbitAverage);
        fprintf(fp, " %e", p.TimeCalcNuclearReactionAverage);
        fprintf(fp, " %e", p.TimeOthersAverage);
        fprintf(fp, "\n");
    }

    static void dumpEachProcess(FILE *fp) {
        WallclockTime & p = getInstance();
        fprintf(fp, "%4d %e\n", PS::Comm::getRank(), p.TimeCalcDensity);
    }

    static void clear() {
        WallclockTime & p = getInstance();
        p.TimeDecomposeDomain      = 0.;
        p.TimeExchangeParticle     = 0.;
        p.TimeCalcDensity          = 0.;
        p.TimeCalcHydro            = 0.;
        p.TimeCalcGravity          = 0.;
        p.TimeReferEquationOfState = 0.;
        p.TimeIntegrateOrbit       = 0.;
        p.TimeCalcNuclearReaction  = 0.;
        p.TimeOthers               = 0.;
    }
};

typedef WallclockTime WT;

namespace PreparedWallclockTime {

    template <class Ttreeforforce>
    void dumpTreeForForce(Ttreeforforce & func,
                          FILE * fp) {
        PS::TimeProfile FuncTime = func.getTimeProfile();
        PS::F64 FuncTimeTotal = FuncTime.getTotalTime();
        fprintf(fp, " %+.2e %+.2e %+.2e %+.2e %+.2e %+.2e %+.2e %+.2e %+.2e %+.2e %+.2e\n",
                FuncTimeTotal,
                FuncTime.make_local_tree,
                FuncTime.make_global_tree,
                FuncTime.calc_force,
                FuncTime.calc_moment_local_tree,
                FuncTime.calc_moment_global_tree,
                FuncTime.make_LET_1st,
                FuncTime.make_LET_2nd,
                FuncTime.exchange_LET_1st,
                FuncTime.exchange_LET_2nd,
                FuncTime.calc_force__core__walk_tree);
        func.clearTimeProfile();
    }

    template <class Tdensity,
              class Thydro,
              class Tgravity>
    void dumpPreparedWallclockTime(PS::F64 time,
                                   FILE * fp,
                                   Tdensity & density,
                                   Thydro & hydro,
                                   Tgravity & gravity) {
        fprintf(fp, "Time: %16.10f\n", time);
        dumpTreeForForce(density, fp);
        dumpTreeForForce(hydro, fp);
        dumpTreeForForce(gravity, fp);
        fprintf(fp, "\n");
    }
};
namespace PWT = PreparedWallclockTime;

