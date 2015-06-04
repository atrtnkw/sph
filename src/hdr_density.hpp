#pragma once

class DensityEPI {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    rs;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
        vel  = sph.vel;
        ksr  = sph.ksr;
        rs   = sph.rs;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
    PS::F64 getRSearch() const {
        return this->rs;
    }
};

class DensityEPJ {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
        vel  = sph.vel;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
};

struct calcDensity {

    void operator () (const DensityEPI *epi,
                      const PS::S32 nip,
                      const DensityEPJ *epj,
                      const PS::S32 njp,
                      Density *density) {
        
        for(PS::S32 i = 0; i < nip; i++) {
            PS::S32    id_i = epi[i].id;
            PS::F64vec x_i  = epi[i].pos;
            PS::F64vec v_i  = epi[i].vel;
            PS::F64    h_i  = epi[i].ksr;

            for(PS::S32 repeat = 0; repeat < 3; repeat++) {
                PS::F64 hi_i  = 1.d / h_i;
                PS::F64 hi3_i = SPH::calcVolumeInverse(hi_i);
                PS::F64 hi4_i = hi_i * hi3_i;
                PS::F64 rh_i  = 0.;
                PS::S32 nj_i  = 0;
                for(PS::S32 j = 0; j < njp; j++) {
                    PS::F64    m_j = epj[j].mass;
                    PS::F64vec x_j = epj[j].pos;
                    
                    PS::F64vec dx_ij = x_i - x_j;
                    PS::F64    r2_ij = dx_ij * dx_ij;
                    PS::F64    r1_ij = sqrt(r2_ij);
                    PS::F64    q_i   = r1_ij * hi_i;
                    
                    PS::F64 kw0 = KernelSph::kernel0th(q_i);
                    
                    PS::F64 rhj =   m_j * hi3_i * kw0;
                    
                    rh_i   += rhj;
                    nj_i   += (q_i < 1.d) ? 1 : 0;
                }
                density[i].dens = rh_i;
                density[i].np   = nj_i;

//                h_i = SPH::eta * KernelSph::ksrh * SPH::calcPowerOfDimInverse(epi[i].mass, rh_i);
                h_i = KernelSph::eta * KernelSph::ksrh
                    * SPH::calcPowerOfDimInverse(epi[i].mass, rh_i);
                density[i].ksr = h_i;
                density[i].itr = (h_i > epi[i].rs) ? true : false;
            }

            PS::F64    hi_i   = 1.d / h_i;
            PS::F64    hi4_i  = hi_i * SPH::calcVolumeInverse(hi_i);
            PS::F64    gh_i  = 0.;
            PS::F64    divv_i = 0.;
            PS::F64vec rotv_i = 0.;            
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S32    id_j = epj[j].id;
                PS::F64    m_j  = epj[j].mass;
                PS::F64vec x_j  = epj[j].pos;
                PS::F64vec v_j  = epj[j].vel;

                PS::F64vec dx_ij = x_i - x_j;
                PS::F64vec dv_ij = v_i - v_j;
                PS::F64    r2_ij = dx_ij * dx_ij;
                PS::F64    r1_ij = sqrt(r2_ij);
                PS::F64    ri_ij = (id_i != id_j) ? 1. / r1_ij : 0.;
                PS::F64    q_i   = r1_ij * hi_i;

                PS::F64 kw0 = KernelSph::kernel0th(q_i);
                PS::F64 kw1 = KernelSph::kernel1st(q_i);
                
                PS::F64    ghj   = - m_j * hi4_i * (KernelSph::dim * kw0 + q_i * kw1);
                PS::F64vec dw_ij = (m_j * hi4_i * kw1 * ri_ij) * dx_ij;

                gh_i   += ghj;
                divv_i += dv_ij * dw_ij;
                rotv_i += dv_ij ^ dw_ij;
            }
            PS::F64 rhi_i = 1. / density[i].dens;
            density[i].grdh = 1.d / (1.d + h_i * rhi_i * gh_i / KernelSph::dim);
//            density[i].grdh = 1.d;
            PS::F64 grd_i = density[i].grdh;
            density[i].rotv = sqrt(rotv_i * rotv_i) * rhi_i * grd_i;
            density[i].divv = divv_i * rhi_i * grd_i;
        }        

    }
};
