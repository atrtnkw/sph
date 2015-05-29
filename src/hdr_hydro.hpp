#pragma once

class DerivativeEPI {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    pres;
    PS::F64    rho;
    PS::F64    bswt;
    PS::F64    vsnd;
    PS::F64    grdh;
    PS::F64    alph;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        pos  = sph.pos;
        vel  = sph.vel;
        ksr  = sph.ksr;
        pres = sph.pres;
        rho  = sph.dens;
        bswt = sph.bswt;
        vsnd = sph.vsnd;
        grdh = sph.grdh;
        alph = sph.alph;
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

class DerivativeEPJ {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    pres;
    PS::F64    rho;
    PS::F64    bswt;
    PS::F64    vsnd;
    PS::F64    grdh;
    PS::F64    alph;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
        vel  = sph.vel;
        ksr  = sph.ksr;
        pres = sph.pres;
        rho  = sph.dens;
        bswt = sph.bswt;
        vsnd = sph.vsnd;
        grdh = sph.grdh;
        alph = sph.alph;
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

struct calcDerivative {

    void operator () (const DerivativeEPI *epi,
                      const PS::S32 nip,
                      const DerivativeEPJ *epj,
                      const PS::S32 njp,
                      Derivative *derivative) {
        for(PS::S32 i = 0; i < nip; i++) {
            PS::S32    id_i    = epi[i].id;
            PS::F64vec x_i     = epi[i].pos;
            PS::F64vec v_i     = epi[i].vel;
            PS::F64    hi_i    = 1. / epi[i].ksr;
            PS::F64    hi4_i   = hi_i * SPH::calcVolumeInverse(hi_i);
            PS::F64    rh_i    = epi[i].rho;
            PS::F64    rhi_i   = 1. / rh_i;
            PS::F64    prhi2_i = epi[i].pres * rhi_i * rhi_i;
            PS::F64    bswt_i  = epi[i].bswt;
            PS::F64    cs_i    = epi[i].vsnd;
            PS::F64    grdh_i  = epi[i].grdh;
            PS::F64    alph_i  = epi[i].alph;
            PS::F64vec acc_i   = 0.;
            PS::F64    ene_i   = 0.;            
            PS::F64    vsmx_i  = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S32    id_j    = epj[j].id;
                PS::F64    m_j     = epj[j].mass;
                PS::F64vec x_j     = epj[j].pos;
                PS::F64vec v_j     = epj[j].vel;
                PS::F64    hi_j    = 1. / epj[j].ksr;
                PS::F64    hi4_j   = hi_j * SPH::calcVolumeInverse(hi_j);
                PS::F64    rh_j    = epj[j].rho;
                PS::F64    rhi_j   = 1. / rh_j;
                PS::F64    prhi2_j = epj[j].pres * rhi_j * rhi_j;
                PS::F64    bswt_j  = epj[j].bswt;
                PS::F64    cs_j    = epj[j].vsnd;
                PS::F64    grdh_j  = epj[j].grdh;
                PS::F64    alph_j  = epj[j].alph;

                PS::F64vec dx_ij = x_i - x_j;
                PS::F64vec dv_ij = v_i - v_j;
                PS::F64    r2_ij = dx_ij * dx_ij;
                PS::F64    r1_ij = sqrt(r2_ij);
                PS::F64    ri_ij = (id_i != id_j) ? 1. / r1_ij : 0.;
                PS::F64    xv_ij = dx_ij * dv_ij;
                PS::F64    q_i   = r1_ij * hi_i;
                PS::F64    q_j   = r1_ij * hi_j;

                PS::F64vec dw_ij  = m_j * dx_ij * ri_ij * 0.5d
                    * (hi4_i * KernelSph::kernelWendlandC2First(q_i)
                       + hi4_j * KernelSph::kernelWendlandC2First(q_j));

                xv_ij = (xv_ij < 0.) ? xv_ij : 0.;
                PS::F64    mu_ij  = xv_ij * ri_ij;
                PS::F64    vs_ij  = cs_i + cs_j - 3. * mu_ij;
                PS::F64    rhi_ij = 2. / (rh_i + rh_j);
                PS::F64    alpha  = 0.5 * (alph_i + alph_j);
                PS::F64    pi_ij  = - 0.5 * vs_ij * mu_ij * rhi_ij * alpha;

                PS::F64    f_ij = 0.5 * (bswt_i + bswt_j);

                acc_i -=         dw_ij * (grdh_i * prhi2_i + grdh_j * prhi2_j + f_ij * pi_ij);
                ene_i += dv_ij * dw_ij * (grdh_i * prhi2_i + 0.5 * f_ij * pi_ij);
                vsmx_i = (vs_ij > vsmx_i) ? vs_ij : vsmx_i;
            }
            derivative[i].acc  = acc_i;
            derivative[i].udot = ene_i;
            derivative[i].vsmx = vsmx_i;
        }

    }
};
