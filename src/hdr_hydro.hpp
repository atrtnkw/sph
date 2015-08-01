#pragma once

class DerivativeEPI {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    pres;
    PS::F64    rho;
    PS::F64    bswt;
    PS::F64    vsnd;
    PS::F64    alph;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        pos  = sph.pos;
        vel  = sph.vel;
        ksr  = sph.ksr;
        hinv = 1.d / sph.ksr;
        pres = sph.pres / (sph.dens * sph.dens) * sph.grdh;
        rho  = sph.dens;
        bswt = sph.bswt;
        vsnd = sph.vsnd;
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
    PS::F64    hinv;
    PS::F64    pres;
    PS::F64    rho;
    PS::F64    bswt;
    PS::F64    vsnd;
    PS::F64    alph;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
        vel  = sph.vel;
        ksr  = sph.ksr;
        hinv = 1.d / sph.ksr;
        pres = sph.pres / (sph.dens * sph.dens) * sph.grdh;
        rho  = sph.dens;
        bswt = sph.bswt;
        vsnd = sph.vsnd;
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

struct calcDerivativeBasic {

    void operator () (const DerivativeEPI *epi,
                      const PS::S32 nip,
                      const DerivativeEPJ *epj,
                      const PS::S32 njp,
                      Derivative *derivative) {
        for(PS::S32 i = 0; i < nip; i++) {
            PS::S32    id_i    = epi[i].id;
            PS::F64vec x_i     = epi[i].pos;
            PS::F64vec v_i     = epi[i].vel;
            PS::F64    hi_i    = epi[i].hinv;
            PS::F64    hi4_i   = hi_i * SPH::calcVolumeInverse(hi_i);
            PS::F64    rh_i    = epi[i].rho;
            PS::F64    prhi2_i = epi[i].pres;
            PS::F64    bswt_i  = epi[i].bswt;
            PS::F64    cs_i    = epi[i].vsnd;
            PS::F64    alph_i  = epi[i].alph;
            PS::F64vec acc_i   = 0.d;
            PS::F64    ene_i   = 0.d;            
            PS::F64    vsmx_i  = 0.d;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S32    id_j    = epj[j].id;
                PS::F64    m_j     = epj[j].mass;
                PS::F64vec x_j     = epj[j].pos;
                PS::F64vec v_j     = epj[j].vel;
                PS::F64    hi_j    = epj[j].hinv;
                PS::F64    hi4_j   = hi_j * SPH::calcVolumeInverse(hi_j);
                PS::F64    rh_j    = epj[j].rho;
                PS::F64    prhi2_j = epj[j].pres;
                PS::F64    bswt_j  = epj[j].bswt;
                PS::F64    cs_j    = epj[j].vsnd;
                PS::F64    alph_j  = epj[j].alph;

                PS::F64vec dx_ij = x_i - x_j;
                PS::F64vec dv_ij = v_i - v_j;
                PS::F64    r2_ij = dx_ij * dx_ij;
                PS::F64    r1_ij = sqrt(r2_ij);
                PS::F64    ri_ij = (id_i != id_j) ? 1.d / r1_ij : 0.d;
                PS::F64    xv_ij = dx_ij * dv_ij;
                PS::F64    q_i   = r1_ij * hi_i;
                PS::F64    q_j   = r1_ij * hi_j;

                PS::F64vec dw_ij  = m_j * dx_ij * ri_ij * 0.5d
                    * (hi4_i * KernelSph::kernel1st(q_i) + hi4_j * KernelSph::kernel1st(q_j));

                PS::F64    mu_ij  = (xv_ij < 0.d) ? (xv_ij * ri_ij) : 0.d;
                PS::F64    vs_ij  = cs_i + cs_j - 3.d * mu_ij;
                PS::F64    pi_ij  = - 0.5d * vs_ij * mu_ij * (alph_i + alph_j) / (rh_i + rh_j);
                PS::F64    f_ij   = 0.5d * (bswt_i + bswt_j);
                PS::F64    vis_ij = f_ij * pi_ij;

                acc_i -=         dw_ij * (prhi2_i + prhi2_j + vis_ij);
                ene_i += dv_ij * dw_ij * (prhi2_i + 0.5d * vis_ij);
                vsmx_i = (vs_ij > vsmx_i) ? vs_ij : vsmx_i;
            }
            derivative[i].acc  = acc_i;
            derivative[i].udot = ene_i;
            derivative[i].vsmx = vsmx_i;
        }

    }
};

struct calcDerivativeX86 {

    void operator () (const DerivativeEPI *epi,
                      const PS::S32 nip,
                      const DerivativeEPJ *epj,
                      const PS::S32 njp,
                      Derivative *derivative) {
#if 0
        PS::S32 nvector = v4df::getVectorLength();
        
        for(PS::S32 i = 0; i < nip; i += nvector) {
            PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;
            for(PS::S32 ii = 0; ii < nii; ii++) {
                PS::S32    id_i    = epi[i+ii].id;
                PS::F64vec x_i     = epi[i+ii].pos;
                PS::F64vec v_i     = epi[i+ii].vel;
                PS::F64    hi_i    = epi[i+ii].hinv;
                PS::F64    hi4_i   = hi_i * SPH::calcVolumeInverse(hi_i);
                PS::F64    rh_i    = epi[i+ii].rho;
                PS::F64    prhi2_i = epi[i+ii].pres;
                PS::F64    bswt_i  = epi[i+ii].bswt;
                PS::F64    cs_i    = epi[i+ii].vsnd;
                PS::F64    alph_i  = epi[i+ii].alph;
                PS::F64vec acc_i   = 0.d;
                PS::F64    ene_i   = 0.d;            
                PS::F64    vsmx_i  = 0.d;
                for(PS::S32 j = 0; j < njp; j++) {
                    PS::S32    id_j    = epj[j].id;
                    PS::F64    m_j     = epj[j].mass;
                    PS::F64vec x_j     = epj[j].pos;
                    PS::F64vec v_j     = epj[j].vel;
                    PS::F64    hi_j    = epj[j].hinv;
                    PS::F64    hi4_j   = hi_j * SPH::calcVolumeInverse(hi_j);
                    PS::F64    rh_j    = epj[j].rho;
                    PS::F64    prhi2_j = epj[j].pres;
                    PS::F64    bswt_j  = epj[j].bswt;
                    PS::F64    cs_j    = epj[j].vsnd;
                    PS::F64    alph_j  = epj[j].alph;
                    
                    PS::F64vec dx_ij = x_i - x_j;
                    PS::F64vec dv_ij = v_i - v_j;
                    PS::F64    r2_ij = dx_ij * dx_ij;
                    PS::F64    r1_ij = sqrt(r2_ij);
                    PS::F64    ri_ij = (id_i != id_j) ? 1.d / r1_ij : 0.d;
                    PS::F64    xv_ij = dx_ij * dv_ij;
                    PS::F64    q_i   = r1_ij * hi_i;
                    PS::F64    q_j   = r1_ij * hi_j;
                    
                    PS::F64vec dw_ij  = m_j * dx_ij * ri_ij * 0.5d
                        * (hi4_i * KernelSph::kernel1st(q_i) + hi4_j * KernelSph::kernel1st(q_j));
                    
                    PS::F64    mu_ij  = (xv_ij < 0.d) ? (xv_ij * ri_ij) : 0.d;
                    PS::F64    vs_ij  = cs_i + cs_j - 3.d * mu_ij;
                    PS::F64    pi_ij  = - 0.5d * vs_ij * mu_ij * (alph_i + alph_j) / (rh_i + rh_j);
                    PS::F64    f_ij   = 0.5d * (bswt_i + bswt_j);
                    PS::F64    vis_ij = f_ij * pi_ij;
                    
                    acc_i -=         dw_ij * (prhi2_i + prhi2_j + vis_ij);
                    ene_i += dv_ij * dw_ij * (prhi2_i + 0.5d * vis_ij);
                    vsmx_i = (vs_ij > vsmx_i) ? vs_ij : vsmx_i;
                }
                derivative[i+ii].acc  = acc_i;
                derivative[i+ii].udot = ene_i;
                derivative[i+ii].vsmx = vsmx_i;
            }
        }
#else
        PS::S32 nvector = v4df::getVectorLength();
        
        for(PS::S32 i = 0; i < nip; i += nvector) {
            v4df id_i(epi[i].id, epi[i+1].id, epi[i+2].id, epi[i+3].id);
            v4df px_i(epi[i].pos[0], epi[i+1].pos[0], epi[i+2].pos[0], epi[i+3].pos[0]);
            v4df py_i(epi[i].pos[1], epi[i+1].pos[1], epi[i+2].pos[1], epi[i+3].pos[1]);
            v4df pz_i(epi[i].pos[2], epi[i+1].pos[2], epi[i+2].pos[2], epi[i+3].pos[2]);
            v4df vx_i(epi[i].vel[0], epi[i+1].vel[0], epi[i+2].vel[0], epi[i+3].vel[0]);
            v4df vy_i(epi[i].vel[1], epi[i+1].vel[1], epi[i+2].vel[1], epi[i+3].vel[1]);
            v4df vz_i(epi[i].vel[2], epi[i+1].vel[2], epi[i+2].vel[2], epi[i+3].vel[2]);
            v4df hi_i(epi[i].hinv, epi[i+1].hinv, epi[i+2].hinv, epi[i+3].hinv);
            v4df rh_i(epi[i].rho, epi[i+1].rho, epi[i+2].rho, epi[i+3].rho);
            v4df prhi2_i(epi[i].pres, epi[i+1].pres, epi[i+2].pres, epi[i+3].pres);
            v4df bswt_i(epi[i].bswt, epi[i+1].bswt, epi[i+2].bswt, epi[i+3].bswt);
            v4df cs_i(epi[i].vsnd, epi[i+1].vsnd, epi[i+2].vsnd, epi[i+3].vsnd);
            v4df alph_i(epi[i].alph, epi[i+1].alph, epi[i+2].alph, epi[i+3].alph);
            v4df hi4_i = hi_i * SPH::calcVolumeInverse(hi_i);
            v4df accx_i(0.d);
            v4df accy_i(0.d);
            v4df accz_i(0.d);
            v4df ene_i(0.d);
            v4df vsmx_i(0.d);

            for(PS::S32 j = 0; j < njp; j++) {
                v4df id_j(epj[j].id);
                v4df m_j(epj[j].mass);
                v4df px_j(epj[j].pos[0]);
                v4df py_j(epj[j].pos[1]);
                v4df pz_j(epj[j].pos[2]);
                v4df vx_j(epj[j].vel[0]);
                v4df vy_j(epj[j].vel[1]);
                v4df vz_j(epj[j].vel[2]);
                v4df hi_j(epj[j].hinv);
                v4df hi4_j = hi_j * SPH::calcVolumeInverse(hi_j);
                v4df rh_j(epj[j].rho);
                v4df prhi2_j(epj[j].pres);
                v4df bswt_j(epj[j].bswt);
                v4df cs_j(epj[j].vsnd);
                v4df alph_j(epj[j].alph);

                v4df dpx_ij = px_i - px_j;
                v4df dpy_ij = py_i - py_j;
                v4df dpz_ij = pz_i - pz_j;
                v4df dvx_ij = vx_i - vx_j;
                v4df dvy_ij = vy_i - vy_j;
                v4df dvz_ij = vz_i - vz_j;

                v4df r2_ij = dpx_ij * dpx_ij;
                r2_ij = v4df::madd(r2_ij, dpy_ij, dpy_ij);
                r2_ij = v4df::madd(r2_ij, dpz_ij, dpz_ij);
                v4df r1_ij = v4df::sqrt(r2_ij);
                v4df ri_ij = v4df(1.d) / r1_ij;
                ri_ij = ((id_i != id_j) & ri_ij);
                v4df xv_ij = dpx_ij * dvx_ij;
                xv_ij = v4df::madd(xv_ij, dpy_ij, dvy_ij);
                xv_ij = v4df::madd(xv_ij, dpz_ij, dvz_ij);
                v4df q_i = r1_ij * hi_i;
                v4df q_j = r1_ij * hi_j;

                v4df dw_ij = v4df(0.5d) * m_j * ri_ij
                    * (hi4_i * KernelSph::kernel1stX(q_i) + hi4_j * KernelSph::kernel1stX(q_j));

                v4df mu_ij = xv_ij * ri_ij;
                mu_ij = ((xv_ij < v4df(0.d)) & mu_ij);
                v4df vs_ij  = cs_i + cs_j - v4df(3.d) * mu_ij;
                v4df pi_ij  = v4df(-0.5d) * vs_ij * mu_ij * (alph_i + alph_j) / (rh_i + rh_j);
                v4df f_ij   = v4df(0.5d) * (bswt_i + bswt_j);
                v4df vis_ij = f_ij * pi_ij;

                v4df da_ij = dw_ij * (prhi2_i + prhi2_j + vis_ij);
                accx_i = v4df::nmadd(accx_i, da_ij, dpx_ij);
                accy_i = v4df::nmadd(accy_i, da_ij, dpy_ij);
                accz_i = v4df::nmadd(accz_i, da_ij, dpz_ij);

                v4df de_ij = v4df::madd(prhi2_i, v4df(0.5d), vis_ij);
                de_ij *= dw_ij;
                ene_i = v4df::madd(ene_i, xv_ij, de_ij);

                vsmx_i = v4df::max(vsmx_i, vs_ij);
            }

            PS::F64 buf_ax[nvector];
            PS::F64 buf_ay[nvector];
            PS::F64 buf_az[nvector];
            PS::F64 buf_eg[nvector];
            PS::F64 buf_vs[nvector];            
            accx_i.store(buf_ax);
            accy_i.store(buf_ay);
            accz_i.store(buf_az);
            ene_i.store(buf_eg);
            vsmx_i.store(buf_vs);

            PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;
            for(PS::S32 ii = 0; ii < nvector; ii++) {
                derivative[i+ii].acc[0] = buf_ax[ii];
                derivative[i+ii].acc[1] = buf_ay[ii];
                derivative[i+ii].acc[2] = buf_az[ii];
                derivative[i+ii].udot   = buf_eg[ii];
                derivative[i+ii].vsmx   = buf_vs[ii];
            }

        }
#endif

    }
};

#ifdef ENABLE_SIMDX86
typedef calcDerivativeX86 calcDerivative;
#else
typedef calcDerivativeBasic calcDerivative;
#endif
