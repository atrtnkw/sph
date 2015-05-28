#pragma once

namespace DampVelocity {    

    const PS::F64 DensityTarget    = 1.0;
    const PS::F64 DensityTolerance = 0.01;

    template <class Tptcl>
    void dampVelocity(Tptcl & system,
                      PS::F64 dt) {
        PS::S32 nloc = system.getNumberOfParticleLocal();
        for(PS::S32 i = 0 ; i < nloc; i++)
            system[i].dampVelocity(dt);    
        
        return;
    }

    template <class Tptcl>
    bool stopDamping(Tptcl & system) {
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

        PS::F64 rho1_var = sqrt(rho2_ave - rho1_ave * rho1_ave);
        if(rho1_var / DensityTarget < DensityTolerance) {
            return true;
        } else {
            return false;
        }
    }
    
}
