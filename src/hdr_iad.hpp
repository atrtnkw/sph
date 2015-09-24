#pragma once

class DerivativeIntegralEPI {
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
    PS::F64mat ctau;
    PS::F64    eta;
    void copyFromFP(const SPH & sph){ 
        id    = sph.id;
        pos   = sph.pos;
        vel   = sph.vel;
        ksr   = sph.ksr;
        hinv  = 1.d / sph.ksr;
        pres  = sph.pres / (sph.dens * sph.dens) * sph.grdh;
        rho   = sph.dens;
        bswt  = 0.5 * sph.bswt;
        vsnd  = sph.vsnd;
        alph  = sph.alph;
        ctau  = sph.ctau;
        eta   = sph.eta;
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

class DerivativeIntegralEPJ {
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
    PS::F64mat ctau;
    PS::F64    eta;
    void copyFromFP(const SPH & sph){ 
        id    = sph.id;
        mass  = sph.mass;
        pos   = sph.pos;
        vel   = sph.vel;
        ksr   = sph.ksr;
        hinv  = 1.d / sph.ksr;
        pres  = sph.pres / (sph.dens * sph.dens) * sph.grdh;
        rho   = sph.dens;
        bswt  = 0.5 * sph.bswt;
        vsnd  = sph.vsnd;
        alph  = sph.alph;
        ctau  = sph.ctau;
        eta   = sph.eta;
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

#ifdef TIMETEST
static PS::F64 tcalchydr = 0.d;
static PS::S32 ncalchydr = 0;
static PS::S64 nintrhydr = 0;
#endif

struct calcDerivativeIntegralBasic {

    void operator () (const DerivativeIntegralEPI *epi,
                      const PS::S32 nip,
                      const DerivativeIntegralEPJ *epj,
                      const PS::S32 njp,
                      Derivative *derivative) {
        for(PS::S32 i = 0; i < nip; i++) {
            PS::S32    id_i    = epi[i].id;
            PS::F64vec x_i     = epi[i].pos;
            PS::F64vec v_i     = epi[i].vel;
            PS::F64    hi_i    = epi[i].hinv;
            PS::F64    hi3_i   = SPH::calcVolumeInverse(hi_i);
            PS::F64    rh_i    = epi[i].rho;
            PS::F64    prhi2_i = epi[i].pres;
            PS::F64    bswt_i  = epi[i].bswt;
            PS::F64    cs_i    = epi[i].vsnd;
            PS::F64    alph_i  = epi[i].alph;
            PS::F64mat ctau_i  = epi[i].ctau;
            PS::F64    eta_i   = epi[i].eta;
            PS::F64vec acc_i   = 0.d;
            PS::F64    ene_i   = 0.d;            
            PS::F64    vsmx_i  = 0.d;
            PS::F64vec g1_i    = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S32    id_j    = epj[j].id;
                PS::F64    m_j     = epj[j].mass;
                PS::F64vec x_j     = epj[j].pos;
                PS::F64vec v_j     = epj[j].vel;
                PS::F64    hi_j    = epj[j].hinv;
                PS::F64    hi3_j   = SPH::calcVolumeInverse(hi_j);
                PS::F64    rh_j    = epj[j].rho;
                PS::F64    prhi2_j = epj[j].pres;
                PS::F64    bswt_j  = epj[j].bswt;
                PS::F64    cs_j    = epj[j].vsnd;
                PS::F64    alph_j  = epj[j].alph;
                PS::F64mat ctau_j  = epj[j].ctau;
                PS::F64    eta_j   = epj[j].eta;

                PS::F64vec dx_ij = x_i - x_j;
                PS::F64vec dv_ij = v_i - v_j;
                PS::F64    r2_ij = dx_ij * dx_ij;
                PS::F64    r1_ij = sqrt(r2_ij);
                PS::F64    ri_ij = (id_j != id_i) ? 1. / r1_ij : 0.;
                PS::F64    xv_ij = dx_ij * dv_ij;
                PS::F64    q_i   = r1_ij * hi_i;
                PS::F64    q_j   = r1_ij * hi_j;
                
                PS::F64vec la_i = 0.;
                PS::F64vec la_j = 0.;
                la_i[0] = ctau_i.xx * dx_ij.x + ctau_i.xy * dx_ij.y + ctau_i.xz * dx_ij.z;
                la_i[1] = ctau_i.xy * dx_ij.x + ctau_i.yy * dx_ij.y + ctau_i.yz * dx_ij.z;
                la_i[2] = ctau_i.xz * dx_ij.x + ctau_i.yz * dx_ij.y + ctau_i.zz * dx_ij.z;
                la_j[0] = ctau_j.xx * dx_ij.x + ctau_j.xy * dx_ij.y + ctau_j.xz * dx_ij.z;
                la_j[1] = ctau_j.xy * dx_ij.x + ctau_j.yy * dx_ij.y + ctau_j.yz * dx_ij.z;
                la_j[2] = ctau_j.xz * dx_ij.x + ctau_j.yz * dx_ij.y + ctau_j.zz * dx_ij.z;
                la_i *= - hi3_i * KernelSph::kernel0th(q_i);
                la_j *=   hi3_j * KernelSph::kernel0th(q_j);
                PS::F64vec la_ij = 0.5 * (la_i - la_j);

                PS::F64 w_ij    = xv_ij * ri_ij;
                PS::F64 w0_ij   = std::min(w_ij, 0.);
                PS::F64 vs_ij   = cs_i + cs_j - 3.d * w0_ij;
                PS::F64 rhi_ij  = 1. / (rh_i + rh_j);
                PS::F64 alph_ij = alph_i + alph_j;
                PS::F64 f_ij    = bswt_i + bswt_j;
                PS::F64 vis0_ij = f_ij * alph_ij * vs_ij;
                PS::F64 vis_ij  = -0.5 * vis0_ij * w0_ij * rhi_ij;

                acc_i -=  m_j * (prhi2_i * la_i - prhi2_j * la_j + vis_ij * la_ij);
                ene_i += (m_j * (prhi2_i * la_i + 0.5 * vis_ij * la_ij)) * dv_ij;
                vsmx_i = (vs_ij > vsmx_i) ? vs_ij : vsmx_i;                                

#ifdef SYMMETRIZED_GRAVITY
                PS::F64vec dg_ij = dx_ij * (m_j * ri_ij * (eta_i * dw_i + eta_j * dw_j));
                g1_i += dg_ij;
#endif
            }
            derivative[i].acc   =       acc_i;
            derivative[i].accg  = 0.5 * g1_i;
            derivative[i].udot  =       ene_i;
            derivative[i].vsmx  = vsmx_i;
        }

    }
};

/*
struct calcDerivativeStandardX86 {

    void operator () (const DerivativeStandardEPI *epi,
                      const PS::S32 nip,
                      const DerivativeStandardEPJ *epj,
                      const PS::S32 njp,
                      Derivative *derivative) {
        v4df (*rcp)(v4df)   = v4df::rcp_4th;
        v4df (*rsqrt)(v4df) = v4df::rsqrt_4th;

        PS::S32 nvector = v4df::getVectorLength();

        for(PS::S32 i = 0; i < nip; i += nvector) {
            v4df id_i(epi[i].id, epi[i+1].id, epi[i+2].id, epi[i+3].id);
            v4df px_i(epi[i].pos[0], epi[i+1].pos[0], epi[i+2].pos[0], epi[i+3].pos[0]);
            v4df py_i(epi[i].pos[1], epi[i+1].pos[1], epi[i+2].pos[1], epi[i+3].pos[1]);
            v4df pz_i(epi[i].pos[2], epi[i+1].pos[2], epi[i+2].pos[2], epi[i+3].pos[2]);
            v4df vx_i(epi[i].vel[0], epi[i+1].vel[0], epi[i+2].vel[0], epi[i+3].vel[0]);
            v4df vy_i(epi[i].vel[1], epi[i+1].vel[1], epi[i+2].vel[1], epi[i+3].vel[1]);
            v4df vz_i(epi[i].vel[2], epi[i+1].vel[2], epi[i+2].vel[2], epi[i+3].vel[2]);
            v4df u_i(epi[i].uene, epi[i+1].uene, epi[i+2].uene, epi[i+3].uene);
            v4df hi_i(epi[i].hinv, epi[i+1].hinv, epi[i+2].hinv, epi[i+3].hinv);
            v4df rh_i(epi[i].rho, epi[i+1].rho, epi[i+2].rho, epi[i+3].rho);
            v4df pres_i(epi[i].pres0, epi[i+1].pres0, epi[i+2].pres0, epi[i+3].pres0);
            v4df prhi2_i(epi[i].pres, epi[i+1].pres, epi[i+2].pres, epi[i+3].pres);
            v4df bswt_i(epi[i].bswt, epi[i+1].bswt, epi[i+2].bswt, epi[i+3].bswt);
            v4df cs_i(epi[i].vsnd, epi[i+1].vsnd, epi[i+2].vsnd, epi[i+3].vsnd);
            v4df alph_i(epi[i].alph, epi[i+1].alph, epi[i+2].alph, epi[i+3].alph);
            v4df alphu_i(epi[i].alphu, epi[i+1].alphu, epi[i+2].alphu, epi[i+3].alphu);
            v4df hi4_i = hi_i * SPH::calcVolumeInverse(hi_i);
            v4df eta_i(epi[i].eta, epi[i+1].eta, epi[i+2].eta, epi[i+3].eta);
            v4df accx_i(0.d);
            v4df accy_i(0.d);
            v4df accz_i(0.d);
            v4df ene_i(0.d);
            v4df vsmx_i(0.d);
            v4df g1x_i(0.d);
            v4df g1y_i(0.d);
            v4df g1z_i(0.d);
            v4df diffu_i(0.);

            for(PS::S32 j = 0; j < njp; j++) {
                v4df dpx_ij = px_i - v4df(epj[j].pos[0]);
                v4df dpy_ij = py_i - v4df(epj[j].pos[1]);
                v4df dpz_ij = pz_i - v4df(epj[j].pos[2]);
                v4df dvx_ij = vx_i - v4df(epj[j].vel[0]);
                v4df dvy_ij = vy_i - v4df(epj[j].vel[1]);
                v4df dvz_ij = vz_i - v4df(epj[j].vel[2]);

                v4df r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v4df ri_ij = rsqrt(r2_ij);
                ri_ij = ((id_i != v4df(epj[j].id)) & ri_ij);
                v4df r1_ij = r2_ij * ri_ij;
                v4df xv_ij = dpx_ij * dvx_ij + dpy_ij * dvy_ij + dpz_ij * dvz_ij;

                v4df hi_j(epj[j].hinv);
                v4df hi4_j = hi_j * SPH::calcVolumeInverse(hi_j);
                v4df q_i = r1_ij * hi_i;
                v4df q_j = r1_ij * hi_j;

                v4df m_j     = v4df(epj[j].mass);
                v4df dw_i    = hi4_i * KernelSph::kernel1st(q_i);
                v4df dw_j    = hi4_j * KernelSph::kernel1st(q_j);
                v4df mdw_ij  = m_j * (dw_i + dw_j);
                v4df mrdw_ij = ri_ij * mdw_ij;

                v4df w_ij    = xv_ij * ri_ij;
                v4df w0_ij   = ((xv_ij < v4df(0.d)) & w_ij);
                v4df vs_ij   = cs_i + v4df(epj[j].vsnd) - v4df(3.d) * w0_ij;
                v4df rhi_ij  = rcp(rh_i + v4df(epj[j].rho));
                v4df alph_ij = alph_i + v4df(epj[j].alph);
                v4df f_ij    = bswt_i + v4df(epj[j].bswt);
                v4df vis0_ij = f_ij * alph_ij * vs_ij;
                v4df vis_ij  = v4df(-0.5d) * vis0_ij * w0_ij * rhi_ij;

                v4df da_ij = mrdw_ij * (prhi2_i + v4df(epj[j].pres) + vis_ij);
                accx_i -= da_ij * dpx_ij;
                accy_i -= da_ij * dpy_ij;
                accz_i -= da_ij * dpz_ij;
                vsmx_i  = v4df::max(vsmx_i, vs_ij);

#ifdef THERMAL_CONDUCTIVITY                
                v4df vsu2_ij  = v4df::fabs(pres_i - v4df(epj[j].pres0)) * rhi_ij * v4df(2.);
                v4df vsui_ij  = ((vsu2_ij != 0.) & rsqrt(vsu2_ij));
                v4df vsu_ij   = vsu2_ij * vsui_ij;
                v4df alphu_ij = alphu_i + v4df(epj[j].alphu);
                v4df u_ij     = u_i - v4df(epj[j].uene);
                v4df de_ij    = prhi2_i * w_ij - rhi_ij
                    * (v4df(0.25) * vis0_ij * w0_ij * w0_ij - alphu_ij * vsu_ij * u_ij);
                ene_i   += mdw_ij * de_ij;
                diffu_i += u_ij * rhi_ij * mrdw_ij;
#else
                v4df de_ij = prhi2_i + v4df(0.5d) * vis_ij;
                de_ij *= mrdw_ij;
                ene_i += xv_ij * de_ij;                
#endif

#ifdef SYMMETRIZED_GRAVITY
                v4df dg_ij  = m_j * ri_ij * (eta_i * dw_i + v4df(epj[j].eta) * dw_j);
                g1x_i += dg_ij * dpx_ij;
                g1y_i += dg_ij * dpy_ij;
                g1z_i += dg_ij * dpz_ij;
#endif
            }

            accx_i  *= v4df(0.5);
            accy_i  *= v4df(0.5);
            accz_i  *= v4df(0.5);
            ene_i   *= v4df(0.5);
            g1x_i   *= v4df(0.5);
            g1y_i   *= v4df(0.5);
            g1z_i   *= v4df(0.5);
            diffu_i *= v4df(2.);
            diffu_i  = v4df::fabs(diffu_i);

            PS::F64 buf_ax[nvector];
            PS::F64 buf_ay[nvector];
            PS::F64 buf_az[nvector];
            PS::F64 buf_eg[nvector];
            PS::F64 buf_vs[nvector];            
            PS::F64 buf_g1x[nvector];            
            PS::F64 buf_g1y[nvector];            
            PS::F64 buf_g1z[nvector];            
            PS::F64 buf_diffu[nvector];
            accx_i.store(buf_ax);
            accy_i.store(buf_ay);
            accz_i.store(buf_az);
            ene_i.store(buf_eg);
            vsmx_i.store(buf_vs);
            g1x_i.store(buf_g1x);
            g1y_i.store(buf_g1y);
            g1z_i.store(buf_g1z);
            diffu_i.store(buf_diffu);

            PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;
            for(PS::S32 ii = 0; ii < nii; ii++) {
                derivative[i+ii].acc[0]  = buf_ax[ii];
                derivative[i+ii].acc[1]  = buf_ay[ii];
                derivative[i+ii].acc[2]  = buf_az[ii];
                derivative[i+ii].udot    = buf_eg[ii];
                derivative[i+ii].vsmx    = buf_vs[ii];
                derivative[i+ii].accg[0] = buf_g1x[ii];
                derivative[i+ii].accg[1] = buf_g1y[ii];
                derivative[i+ii].accg[2] = buf_g1z[ii];
                derivative[i+ii].diffu   = buf_diffu[ii];
            }

        }

    }
};
*/

#ifdef INTEGRAL_APPROACH_DERIVATIVE
typedef DerivativeIntegralEPI DerivativeEPI;
typedef DerivativeIntegralEPJ DerivativeEPJ;
//#ifdef ENABLE_SIMDX86
//typedef calcDerivativeIntegralX86 calcDerivative;
//#elif defined ENABLE_SIMDX86_SINGLE
//typedef calcDerivativeIntegralSingleX86 calcDerivative;
//#else
typedef calcDerivativeIntegralBasic calcDerivative;
//#endif
#endif
