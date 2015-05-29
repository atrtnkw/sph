//////////////////////////////////////////////////////////////////
// 'MaxOfSmallVariance' contains some bugs.
//////////////////////////////////////////////////////////////////

#pragma once

namespace DampVelocity {    

    const PS::F64 DensityTarget    = 1.0;
    const PS::F64 DensityTolerance = 1e-4;
    const PS::F64 MaxOfSmallVariance = 1e-4;
    const PS::S32 NumberOfSmallVariance = 1000;

    bool first = true;
    FILE *fp;
    PS::F64 rho1_var_prev = 1e30;
    PS::S32 n_var_dev = 0;

    template <class Tptcl>
    void dampVelocity(Tptcl & system,
                      PS::F64 dt) {
        PS::S32 nloc = system.getNumberOfParticleLocal();
        for(PS::S32 i = 0 ; i < nloc; i++)
            system[i].dampVelocity(dt);    
        
        return;
    }

    template <class Tptcl>
    bool stopDamping(PS::F64 time,
                     Tptcl & system) {
        PS::F64 rho1_loc = 0.0;
        PS::F64 rho2_loc = 0.0;

        PS::S32 nloc = system.getNumberOfParticleLocal();
        for(PS::S32 i = 0 ; i < nloc; i++) {
            rho1_loc += system[i].dens;
            rho2_loc += system[i].dens * system[i].dens;
        }        

        PS::S32 nglb     = PS::Comm::getSum(nloc);
        PS::F64 rho1_glb = PS::Comm::getSum(rho1_loc);
        PS::F64 rho2_glb = PS::Comm::getSum(rho2_loc);

        PS::F64 rho1_ave = rho1_glb / (PS::S32)nglb;
        PS::F64 rho2_ave = rho2_glb / (PS::S32)nglb;

        PS::F64 rho1_var = sqrt(rho2_ave - rho1_ave * rho1_ave) / DensityTarget;
        PS::F64 var_dev  = fabs((rho1_var - rho1_var_prev) / rho1_var_prev);

        if(first) {
            fp = fopen("snap/damping.log", "w");
            first = false;
        }        
        fprintf(fp, "time: %.10f %+e %+e\n", time, rho1_var, var_dev);
        fflush(fp);

        if(var_dev < MaxOfSmallVariance) {
            n_var_dev++;
        }

        if(rho1_var < DensityTolerance || n_var_dev > NumberOfSmallVariance) {
            rho1_var_prev = rho1_var;
            return true;
        } else {
            rho1_var_prev = rho1_var;
            return false;
        }
    }
    
}
