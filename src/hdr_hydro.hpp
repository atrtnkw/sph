#pragma once

class DerivativeStandardEPI {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    uene;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    pres0;
    PS::F64    pres;
    PS::F64    rho;
    PS::F64    bswt;
    PS::F64    vsnd;
    PS::F64    alph;
    PS::F64    alphu;
    PS::F64    eta;
    void copyFromFP(const SPH & sph){ 
        id    = sph.id;
        pos   = sph.pos;
        vel   = sph.vel;
        ksr   = sph.ksr;
        hinv  = 1.d / sph.ksr;
        pres0 = sph.pres;
        pres  = sph.pres / (sph.dens * sph.dens) * sph.grdh;
        rho   = sph.dens;
        bswt  = 0.5 * sph.bswt;
        vsnd  = sph.vsnd;
        alph  = sph.alph;
        alphu = sph.alphu;
        eta   = sph.eta;
        PS::F64 umin = CalcEquationOfState::getEnergyMin(sph.dens, sph.abar, sph.zbar);
        uene  = std::max(sph.uene - umin, 0.);
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

class DerivativeStandardEPJ {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    uene;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    pres0;
    PS::F64    pres;
    PS::F64    rho;
    PS::F64    bswt;
    PS::F64    vsnd;
    PS::F64    alph;
    PS::F64    alphu;
    PS::F64    eta;
    void copyFromFP(const SPH & sph){ 
        id    = sph.id;
        mass  = sph.mass;
        pos   = sph.pos;
        vel   = sph.vel;
        //uene  = sph.uene;
        ksr   = sph.ksr;
        hinv  = 1.d / sph.ksr;
        pres0 = sph.pres;
        pres  = sph.pres / (sph.dens * sph.dens) * sph.grdh;
        rho   = sph.dens;
        bswt  = 0.5 * sph.bswt;
        vsnd  = sph.vsnd;
        alph  = sph.alph;
        alphu = sph.alphu;
        eta   = sph.eta;
        PS::F64 umin = CalcEquationOfState::getEnergyMin(sph.dens, sph.abar, sph.zbar);
        uene  = std::max(sph.uene - umin, 0.);
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

struct calcDerivativeStandardBasic {

    void operator () (const DerivativeStandardEPI *epi,
                      const PS::S32 nip,
                      const DerivativeStandardEPJ *epj,
                      const PS::S32 njp,
                      Derivative *derivative) {
        for(PS::S32 i = 0; i < nip; i++) {
            PS::S32    id_i    = epi[i].id;
            PS::F64vec x_i     = epi[i].pos;
            PS::F64vec v_i     = epi[i].vel;
            PS::F64    u_i     = epi[i].uene;
            PS::F64    hi_i    = epi[i].hinv;
            PS::F64    hi4_i   = hi_i * SPH::calcVolumeInverse(hi_i);
            PS::F64    rh_i    = epi[i].rho;
            PS::F64    pres_i  = epi[i].pres0;
            PS::F64    prhi2_i = epi[i].pres;
            PS::F64    bswt_i  = epi[i].bswt;
            PS::F64    cs_i    = epi[i].vsnd;
            PS::F64    alph_i  = epi[i].alph;
            PS::F64    alphu_i = epi[i].alphu;
            PS::F64    eta_i   = epi[i].eta;
            PS::F64vec acc_i   = 0.d;
            PS::F64    ene_i   = 0.d;            
            PS::F64    vsmx_i  = 0.d;
            PS::F64vec g1_i    = 0.;
            PS::F64    diffu_i = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S32    id_j    = epj[j].id;
                PS::F64    m_j     = epj[j].mass;
                PS::F64vec x_j     = epj[j].pos;
                PS::F64vec v_j     = epj[j].vel;
                PS::F64    u_j     = epj[j].uene;
                PS::F64    hi_j    = epj[j].hinv;
                PS::F64    hi4_j   = hi_j * SPH::calcVolumeInverse(hi_j);
                PS::F64    rh_j    = epj[j].rho;
                PS::F64    pres_j  = epj[j].pres0;
                PS::F64    prhi2_j = epj[j].pres;
                PS::F64    bswt_j  = epj[j].bswt;
                PS::F64    cs_j    = epj[j].vsnd;
                PS::F64    alph_j  = epj[j].alph;
                PS::F64    alphu_j = epj[j].alphu;
                PS::F64    eta_j   = epj[j].eta;

                PS::F64vec dx_ij = x_i - x_j;
                PS::F64vec dv_ij = v_i - v_j;
                PS::F64    r2_ij = dx_ij * dx_ij;
                PS::F64    r1_ij = sqrt(r2_ij);
                PS::F64    ri_ij = (id_i != id_j) ? 1.d / r1_ij : 0.d;
                PS::F64    xv_ij = dx_ij * dv_ij;
                PS::F64    q_i   = r1_ij * hi_i;
                PS::F64    q_j   = r1_ij * hi_j;

                PS::F64 dw_i    = hi4_i * KernelSph::kernel1st(q_i);
                PS::F64 dw_j    = hi4_j * KernelSph::kernel1st(q_j);
                PS::F64 mdw_ij  = m_j * (dw_i + dw_j);
                PS::F64 mrdw_ij = ri_ij * mdw_ij;

                PS::F64 w_ij    = xv_ij * ri_ij;
                PS::F64 w0_ij   = std::min(w_ij, 0.);
                PS::F64 vs_ij   = cs_i + cs_j - 3.d * w0_ij;
                PS::F64 rhi_ij  = 1. / (rh_i + rh_j);
                PS::F64 alph_ij = alph_i + alph_j;
                PS::F64 f_ij    = bswt_i + bswt_j;
                PS::F64 vis0_ij = f_ij * alph_ij * vs_ij;
                PS::F64 vis_ij  = -0.5 * vis0_ij * w0_ij * rhi_ij;

                acc_i -= (mrdw_ij * (prhi2_i + prhi2_j + vis_ij)) * dx_ij;
                vsmx_i = (vs_ij > vsmx_i) ? vs_ij : vsmx_i;                                

#ifdef THERMAL_CONDUCTIVITY
                PS::F64 vsu_ij   = sqrt(fabs(pres_i - pres_j) * rhi_ij * 2.);
                PS::F64 alphu_ij = alphu_i + alphu_j;
                PS::F64 u_ij     = u_i - u_j;
                PS::F64 de_ij    = prhi2_i * w_ij - rhi_ij
                    * (0.25 * vis0_ij * w0_ij * w0_ij - alphu_ij * vsu_ij * u_ij);
                ene_i   += mdw_ij * de_ij;
                diffu_i += u_ij * rhi_ij * mrdw_ij;
#else
                PS::F64 de_ij = prhi2_i + 0.5 * vis_ij;
                ene_i += de_ij * mdw_ij * w_ij;
#endif

#ifdef SYMMETRIZED_GRAVITY
                PS::F64vec dg_ij = dx_ij * (m_j * ri_ij * (eta_i * dw_i + eta_j * dw_j));
                g1_i += dg_ij;
#endif
            }
            derivative[i].acc   = 0.5 * acc_i;
            derivative[i].accg  = 0.5 * g1_i;
            derivative[i].udot  = 0.5 * ene_i;
            derivative[i].vsmx  = vsmx_i;
            derivative[i].diffu = 4. * fabs(diffu_i);
        }

    }
};

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
            PS::F64 buf_id[nvector];
            PS::F64 buf_px[nvector];
            PS::F64 buf_py[nvector];
            PS::F64 buf_pz[nvector];
            PS::F64 buf_vx[nvector];
            PS::F64 buf_vy[nvector];
            PS::F64 buf_vz[nvector];
            PS::F64 buf_u[nvector];
            PS::F64 buf_hi[nvector];
            PS::F64 buf_rh[nvector];
            PS::F64 buf_pres[nvector];
            PS::F64 buf_prhi2[nvector];
            PS::F64 buf_bswt[nvector];
            PS::F64 buf_cs[nvector];
            PS::F64 buf_alph[nvector];
            PS::F64 buf_alphu[nvector];
            PS::F64 buf_eta[nvector];
            PS::S32 nii = std::min(nip - i, nvector);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_id[ii]    = epi[i+ii].id;
                buf_px[ii]    = epi[i+ii].pos[0];
                buf_py[ii]    = epi[i+ii].pos[1];
                buf_pz[ii]    = epi[i+ii].pos[2];
                buf_vx[ii]    = epi[i+ii].vel[0];
                buf_vy[ii]    = epi[i+ii].vel[1];
                buf_vz[ii]    = epi[i+ii].vel[2];
                buf_u[ii]     = epi[i+ii].uene;
                buf_hi[ii]    = epi[i+ii].hinv;
                buf_rh[ii]    = epi[i+ii].rho;
                buf_pres[ii]  = epi[i+ii].pres0;
                buf_prhi2[ii] = epi[i+ii].pres;
                buf_bswt[ii]  = epi[i+ii].bswt;
                buf_cs[ii]    = epi[i+ii].vsnd;
                buf_alph[ii]  = epi[i+ii].alph;
                buf_alphu[ii] = epi[i+ii].alphu;
                buf_eta[ii]   = epi[i+ii].eta;
            }
            v4df id_i;
            v4df px_i, py_i, pz_i;
            v4df vx_i, vy_i, vz_i;
            v4df u_i,  hi_i, rh_i;
            v4df pres_i, prhi2_i, bswt_i;
            v4df cs_i, alph_i, alphu_i;
            v4df eta_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            u_i.load(buf_u);
            hi_i.load(buf_hi);
            rh_i.load(buf_rh);
            pres_i.load(buf_pres);
            prhi2_i.load(buf_prhi2);
            bswt_i.load(buf_bswt);
            cs_i.load(buf_cs);
            alph_i.load(buf_alph);
            alphu_i.load(buf_alphu);
            eta_i.load(buf_eta);
            v4df hi4_i = hi_i * SPH::calcVolumeInverse(hi_i);
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

            //PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;
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

#if 0
struct calcDerivativeSingleX86 {

#ifdef SYMMETRIZED_GRAVITY
#error SYMMETRIZED_GRAVITY has been not yet implemented in float precision!
#endif
#ifdef SYMMETRIZED_GRAVITY
#error Not adjusted to extended to IAD!
#endif
    void operator () (const DerivativeEPI *epi,
                      const PS::S32 nip,
                      const DerivativeEPJ *epj,
                      const PS::S32 njp,
                      Derivative *derivative) {
#ifdef TIMETEST
        nintrhydr += nip * njp;
#endif
        v8sf (*rcp)(v8sf)   = v8sf::rcp_1st;
        v8sf (*rsqrt)(v8sf) = v8sf::rsqrt_1st;

        PS::S32 nvector = v8sf::getVectorLength();
        
        for(PS::S32 i = 0; i < nip; i += nvector) {
            v8sf id_i(epi[i].id,   epi[i+1].id, epi[i+2].id, epi[i+3].id,
                      epi[i+4].id, epi[i+5].id, epi[i+6].id, epi[i+7].id);
            v8sf px_i(epi[i].pos[0],   epi[i+1].pos[0], epi[i+2].pos[0], epi[i+3].pos[0],
                      epi[i+4].pos[0], epi[i+5].pos[0], epi[i+6].pos[0], epi[i+7].pos[0]);
            v8sf py_i(epi[i].pos[1],   epi[i+1].pos[1], epi[i+2].pos[1], epi[i+3].pos[1],
                      epi[i+4].pos[1], epi[i+5].pos[1], epi[i+6].pos[1], epi[i+7].pos[1]);
            v8sf pz_i(epi[i].pos[2],   epi[i+1].pos[2], epi[i+2].pos[2], epi[i+3].pos[2],
                      epi[i+4].pos[2], epi[i+5].pos[2], epi[i+6].pos[2], epi[i+7].pos[2]);
            v8sf vx_i(epi[i].vel[0],   epi[i+1].vel[0], epi[i+2].vel[0], epi[i+3].vel[0],
                      epi[i+4].vel[0], epi[i+5].vel[0], epi[i+6].vel[0], epi[i+7].vel[0]);
            v8sf vy_i(epi[i].vel[1],   epi[i+1].vel[1], epi[i+2].vel[1], epi[i+3].vel[1],
                      epi[i+4].vel[1], epi[i+5].vel[1], epi[i+6].vel[1], epi[i+7].vel[1]);
            v8sf vz_i(epi[i].vel[2],   epi[i+1].vel[2], epi[i+2].vel[2], epi[i+3].vel[2],
                      epi[i+4].vel[2], epi[i+5].vel[2], epi[i+6].vel[2], epi[i+7].vel[2]);
            v8sf hi_i(epi[i].hinv,   epi[i+1].hinv, epi[i+2].hinv, epi[i+3].hinv,
                      epi[i+4].hinv, epi[i+5].hinv, epi[i+6].hinv, epi[i+7].hinv);
            v8sf rh_i(epi[i].rho,   epi[i+1].rho, epi[i+2].rho, epi[i+3].rho,
                      epi[i+4].rho, epi[i+5].rho, epi[i+6].rho, epi[i+7].rho);
            v8sf prhi2_i(epi[i].pres,   epi[i+1].pres, epi[i+2].pres, epi[i+3].pres,
                         epi[i+4].pres, epi[i+5].pres, epi[i+6].pres, epi[i+7].pres);
            v8sf bswt_i(epi[i].bswt,   epi[i+1].bswt, epi[i+2].bswt, epi[i+3].bswt,
                        epi[i+4].bswt, epi[i+5].bswt, epi[i+6].bswt, epi[i+7].bswt);
            v8sf cs_i(epi[i].vsnd,   epi[i+1].vsnd, epi[i+2].vsnd, epi[i+3].vsnd,
                      epi[i+4].vsnd, epi[i+5].vsnd, epi[i+6].vsnd, epi[i+7].vsnd);
            v8sf alph_i(epi[i].alph,   epi[i+1].alph, epi[i+2].alph, epi[i+3].alph,
                        epi[i+4].alph, epi[i+5].alph, epi[i+6].alph, epi[i+7].alph);
            v8sf hi4_i = hi_i * SPH::calcVolumeInverse(hi_i);
            v8sf accx_i(0.d);
            v8sf accy_i(0.d);
            v8sf accz_i(0.d);
            v8sf ene_i(0.d);
            v8sf vsmx_i(0.d);

#ifdef TIMETEST
            PS::F64 t1 = getWallclockTime();
#endif

            for(PS::S32 j = 0; j < njp; j++) {
                v8sf dpx_ij = px_i - v8sf(epj[j].pos[0]);
                v8sf dpy_ij = py_i - v8sf(epj[j].pos[1]);
                v8sf dpz_ij = pz_i - v8sf(epj[j].pos[2]);
                v8sf dvx_ij = vx_i - v8sf(epj[j].vel[0]);
                v8sf dvy_ij = vy_i - v8sf(epj[j].vel[1]);
                v8sf dvz_ij = vz_i - v8sf(epj[j].vel[2]);

                v8sf r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v8sf ri_ij = rsqrt(r2_ij);
                ri_ij = ((id_i != v8sf(epj[j].id)) & ri_ij);
                v8sf r1_ij = r2_ij * ri_ij;
                v8sf xv_ij = dpx_ij * dvx_ij + dpy_ij * dvy_ij + dpz_ij * dvz_ij;

                v8sf hi_j(epj[j].hinv);
                v8sf hi4_j = hi_j * SPH::calcVolumeInverse(hi_j);
                v8sf q_i = r1_ij * hi_i;
                v8sf q_j = r1_ij * hi_j;

                v8sf dw_ij = v8sf(0.5d) * v8sf(epj[j].mass) * ri_ij
                    * (hi4_i * KernelSph::kernel1st(q_i) + hi4_j * KernelSph::kernel1st(q_j));

                v8sf mu_ij = xv_ij * ri_ij;
                mu_ij = ((xv_ij < v8sf(0.d)) & mu_ij);
                v8sf vs_ij  = cs_i + v8sf(epj[j].vsnd) - v8sf(3.d) * mu_ij;
                v8sf pi_ij  = v8sf(-0.5d) * vs_ij * mu_ij * (alph_i + v8sf(epj[j].alph))
                    * rcp(rh_i + v8sf(epj[j].rho));
                v8sf f_ij   = v8sf(0.5d) * (bswt_i + v8sf(epj[j].bswt));
                v8sf vis_ij = f_ij * pi_ij;

                v8sf da_ij = dw_ij * (prhi2_i + v8sf(epj[j].pres) + vis_ij);
                accx_i -= da_ij * dpx_ij;
                accy_i -= da_ij * dpy_ij;
                accz_i -= da_ij * dpz_ij;

                v8sf de_ij = prhi2_i + v8sf(0.5d) * vis_ij;
                de_ij *= dw_ij;
                ene_i += xv_ij * de_ij;

                vsmx_i = v8sf::max(vsmx_i, vs_ij);
            }
#ifdef TIMETEST
            tcalchydr += getWallclockTime() - t1;
            ncalchydr++;
#endif

            PS::F32 buf_ax[nvector];
            PS::F32 buf_ay[nvector];
            PS::F32 buf_az[nvector];
            PS::F32 buf_eg[nvector];
            PS::F32 buf_vs[nvector];            
            accx_i.store(buf_ax);
            accy_i.store(buf_ay);
            accz_i.store(buf_az);
            ene_i.store(buf_eg);
            vsmx_i.store(buf_vs);

            PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;
            //for(PS::S32 ii = 0; ii < nvector; ii++) {
            for(PS::S32 ii = 0; ii < nii; ii++) {
                derivative[i+ii].acc[0] = buf_ax[ii];
                derivative[i+ii].acc[1] = buf_ay[ii];
                derivative[i+ii].acc[2] = buf_az[ii];
                derivative[i+ii].udot   = buf_eg[ii];
                derivative[i+ii].vsmx   = buf_vs[ii];
            }

        }

    }
};
#endif

#ifndef INTEGRAL_APPROACH_DERIVATIVE
typedef DerivativeStandardEPI DerivativeEPI;
typedef DerivativeStandardEPJ DerivativeEPJ;
#ifdef ENABLE_SIMDX86
typedef calcDerivativeStandardX86 calcDerivative;
#elif defined ENABLE_SIMDX86_SINGLE
typedef calcDerivativeStandardSingleX86 calcDerivative;
#else
typedef calcDerivativeStandardBasic calcDerivative;
#endif
#endif
