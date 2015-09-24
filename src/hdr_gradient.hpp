#pragma once

class GradientEPI {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64    ksr;
    void copyFromFP(const SPH & sph) {
        id  = sph.id;
        pos = sph.pos;
        ksr = sph.ksr;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
    PS::F64 getRSearch() const {
        return this->ksr;
    }    
};

class GradientEPJ {
public:
    PS::S32    id;
    PS::F64    vol;
    PS::F64vec pos;
    void copyFromFP(const SPH & sph) {
        id  = sph.id;
        vol = sph.mass / sph.dens;
        pos = sph.pos;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
};

struct calcGradientBasic {

    static PS::F64mat invertMatrix(PS::F64mat & tau) {
        PS::F64mat c = 0.;
#ifdef USE_AT1D
        c.xx = 1. / tau.xx;
#elif defined USE_AT2D
        PS::F64 detinv = 1. / (tau.xx * tau.yy - tau.xy * tau.xy);
        c.xx =   tau.yy * detinv;
        c.yy =   tau.xx * detinv;
        c.xy = - tau.xy * detinv;
#else
        PS::F64 detinv = 1. / (tau.xx * tau.yy * tau.zz
                               + tau.xy * tau.xz * tau.yz * 2.
                               - tau.xx * tau.yz * tau.yz
                               - tau.yy * tau.xz * tau.xz
                               - tau.zz * tau.xy * tau.xy);
        c.xx = (tau.yy * tau.zz - tau.yz * tau.yz) * detinv;
        c.yy = (tau.xx * tau.zz - tau.xz * tau.xz) * detinv;
        c.zz = (tau.xx * tau.yy - tau.xy * tau.xy) * detinv;
        c.xy = (tau.xz * tau.yz - tau.xy * tau.zz) * detinv;
        c.xz = (tau.xy * tau.yz - tau.xz * tau.yy) * detinv;
        c.yz = (tau.xz * tau.xy - tau.xx * tau.yz) * detinv;
#endif
        return c;
    }

    void operator () (const GradientEPI *epi,
                      const PS::S32 nip,
                      const GradientEPJ *epj,
                      const PS::S32 njp,
                      Gradient * gradient) {

        for(PS::S32 i = 0; i < nip; i++) {
            PS::S32    id_i  = epi[i].id;
            PS::F64vec x_i   = epi[i].pos;
            PS::F64    hi_i  = 1. / epi[i].ksr;
            PS::F64    hi3_i = SPH::calcVolumeInverse(hi_i);
            PS::F64mat tau_i = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S32    id_j = epj[j].id;
                PS::F64    vl_j = epj[j].vol;
                PS::F64vec x_j  = epj[j].pos;

                PS::F64vec dx_ij = x_i - x_j;
                PS::F64    r2_ij = dx_ij * dx_ij;
                PS::F64    r1_ij = sqrt(r2_ij);
                PS::F64    q_i   = r1_ij * hi_i;
                PS::F64    vw0_i = vl_j * hi3_i * KernelSph::kernel0th(q_i);

                tau_i.xx += dx_ij[0] * dx_ij[0] * vw0_i;
                tau_i.yy += dx_ij[1] * dx_ij[1] * vw0_i;
                tau_i.zz += dx_ij[2] * dx_ij[2] * vw0_i;
                tau_i.xy += dx_ij[0] * dx_ij[1] * vw0_i;
                tau_i.xz += dx_ij[0] * dx_ij[2] * vw0_i;
                tau_i.yz += dx_ij[1] * dx_ij[2] * vw0_i;
            }
            gradient[i].ctau = invertMatrix(tau_i);
        }
    }
};

typedef calcGradientBasic calcGradient;
