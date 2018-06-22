#pragma once

class Hydro {
public:
    PS::F64vec acch;
    PS::F64vec accg;
    PS::F64    udot;
    PS::F64    vsmx;
    PS::F64    diffu;
    void clear() {
        this->acch  = 0.;
        this->accg  = 0.;
        this->udot  = 0.;
        this->vsmx  = 0.;
        this->diffu = 0.;
    }
};

void SPH::copyFromForce(const Hydro & hydro) {
    using namespace CodeUnit;
    this->acch  = hydro.acch;
    this->accg1 = hydro.accg * GravityConstantInThisUnit;
    this->udot  = hydro.udot;
    this->vsmx  = hydro.vsmx;
    this->diffu = hydro.diffu;
}

class HydroEPI {
public:
    PS::S64    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    pres;
    PS::F64    presc;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;
    PS::F64    alphu;
    PS::F64    thrm;
    PS::F64    eta;

    void copyFromFP(const SPH & sph){ 
        this->id    = sph.id;
        this->pos   = sph.pos;
        this->vel   = sph.vel;
        this->ksr   = sph.ksr;
        this->hinv  = 1. / sph.ksr;
        this->pres  = sph.pres;
        this->presc = sph.pres * sph.grdh / (sph.dens * sph.dens);
        this->bswt  = 0.5 * sph.bswt;
        this->dens  = sph.dens;
        this->vsnd  = sph.vsnd;
        this->alph  = sph.alph;
        this->alphu = sph.alphu;        
        //PS::F64 ene = sph.uene - CalcEquationOfState::getEnergyMin(sph.dens, sph.abar, sph.zbar);
        PS::F64 ene = sph.uene - sph.umin;
        this->thrm  = std::max(ene, 0.);
        this->eta   = sph.eta;
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

class HydroEPJ {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    pres;
    PS::F64    presc;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;
    PS::F64    alphu;
    PS::F64    thrm;
    PS::F64    eta;

    void copyFromFP(const SPH & sph){ 
        this->id    = sph.id;
        this->mass  = sph.mass;
        this->pos   = sph.pos;
        this->vel   = sph.vel;
        this->ksr   = sph.ksr;
        this->hinv  = 1. / sph.ksr;
        this->pres  = sph.pres;
        this->presc = sph.pres * sph.grdh / (sph.dens * sph.dens);
        this->bswt  = 0.5 * sph.bswt;
        this->dens  = sph.dens;
        this->vsnd  = sph.vsnd;
        this->alph  = sph.alph;
        this->alphu = sph.alphu;
        //PS::F64 ene = sph.uene - CalcEquationOfState::getEnergyMin(sph.dens, sph.abar, sph.zbar);
        PS::F64 ene = sph.uene - sph.umin;
        this->thrm  = std::max(ene, 0.);
        this->eta   = sph.eta;
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

#ifdef USE_INTRINSICS

struct calcHydro {

    void operator () (const HydroEPI * epi,
                      const PS::S32 nip,
                      const HydroEPJ * epj,
                      const PS::S32 njp,
                      Hydro * hydro) {
        vndf (*rcp)(vndf)   = vndf::rcp_4th;
        vndf (*rsqrt)(vndf) = vndf::rsqrt_4th;

        PS::S32 nvector = vndf::getVectorLength();

        for(PS::S32 i = 0; i < nip; i += nvector) {
            PS::F64 buf_id[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_px[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_py[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_pz[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_vx[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_vy[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_vz[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_hi[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_pp[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_pc[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_bs[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_dn[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_cs[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_al[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_au[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_eg[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_et[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::S32 nii = std::min(nip - i, nvector);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_id[ii] = epi[i+ii].id;
                buf_px[ii] = epi[i+ii].pos[0];
                buf_py[ii] = epi[i+ii].pos[1];
                buf_pz[ii] = epi[i+ii].pos[2];
                buf_vx[ii] = epi[i+ii].vel[0];
                buf_vy[ii] = epi[i+ii].vel[1];
                buf_vz[ii] = epi[i+ii].vel[2];
                buf_hi[ii] = epi[i+ii].hinv;
                buf_pp[ii] = epi[i+ii].pres;
                buf_pc[ii] = epi[i+ii].presc;
                buf_bs[ii] = epi[i+ii].bswt;
                buf_dn[ii] = epi[i+ii].dens;
                buf_cs[ii] = epi[i+ii].vsnd;
                buf_al[ii] = epi[i+ii].alph;
                buf_au[ii] = epi[i+ii].alphu;
                buf_eg[ii] = epi[i+ii].thrm;
                buf_et[ii] = epi[i+ii].eta;
            }
            vndf id_i;
            vndf px_i, py_i, pz_i;
            vndf vx_i, vy_i, vz_i;
            vndf hi_i, pp_i, pc_i;
            vndf bs_i, dn_i, cs_i, al_i, au_i, eg_i;
            vndf et_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            hi_i.load(buf_hi);
            pp_i.load(buf_pp);
            pc_i.load(buf_pc);
            bs_i.load(buf_bs);
            dn_i.load(buf_dn);
            cs_i.load(buf_cs);
            al_i.load(buf_al);
            au_i.load(buf_au);
            eg_i.load(buf_eg);
            et_i.load(buf_et);
            vndf hi4_i = hi_i * ND::calcVolumeInverse(hi_i);
            vndf achx_i(0.);
            vndf achy_i(0.);
            vndf achz_i(0.);
            vndf udot_i(0.);
            vndf vsmx_i(0.);
            vndf difu_i(0.);
            vndf g1x_i(0.);
            vndf g1y_i(0.);
            vndf g1z_i(0.);
            for(PS::S32 j = 0; j < njp; j++) {
                vndf dpx_ij = px_i - vndf(epj[j].pos[0]);
                vndf dpy_ij = py_i - vndf(epj[j].pos[1]);
                vndf dpz_ij = pz_i - vndf(epj[j].pos[2]);
                vndf dvx_ij = vx_i - vndf(epj[j].vel[0]);
                vndf dvy_ij = vy_i - vndf(epj[j].vel[1]);
                vndf dvz_ij = vz_i - vndf(epj[j].vel[2]);

                vndf r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
#ifdef __AVX512F__
                vndf ri_ij = _mm512_mask_mov_pd(vndf(0.),
                    _mm512_cmp_pd_mask(id_i, vndf(epj[j].id), _CMP_NEQ_OQ),
                    rsqrt(r2_ij));
#else
                vndf ri_ij = ((id_i != vndf(epj[j].id)) & rsqrt(r2_ij));
#endif
                vndf r1_ij = r2_ij * ri_ij;
                vndf q_i   = r1_ij * hi_i;
                vndf q_j   = r1_ij * epj[j].hinv;

                vndf hi4_j = vndf(epj[j].hinv) * ND::calcVolumeInverse(epj[j].hinv);
                vndf dw_i  = hi4_i * SK::kernel1st(q_i);
                vndf dw_j  = hi4_j * SK::kernel1st(q_j);
                vndf ka_ij = (dw_i + dw_j) * vndf(epj[j].mass);

                vndf rv_ij  = dpx_ij * dvx_ij + dpy_ij * dvy_ij + dpz_ij * dvz_ij;
                vndf w_ij   = rv_ij * ri_ij;
#ifdef __AVX512F__
                vndf w0_ij  = _mm512_mask_mov_pd(vndf(0.),
                    _mm512_cmp_pd_mask(w_ij, vndf(0.), _CMP_LT_OS),
                    w_ij);
#else
                vndf w0_ij  = ((w_ij < vndf(0.)) & w_ij);
#endif
                vndf vs_ij  = cs_i + vndf(epj[j].vsnd) - vndf(3.) * w0_ij;
                vndf rhi_ij = rcp(dn_i + vndf(epj[j].dens));
                vndf av0_ij = (bs_i + vndf(epj[j].bswt)) * (al_i + vndf(epj[j].alph))
                    * vs_ij * w0_ij;
                                
                vsmx_i  = vndf::max(vsmx_i, vs_ij);

                vndf ta_ij = (pc_i + vndf(epj[j].presc) - vndf(0.5) * av0_ij * rhi_ij)
                    * ka_ij * ri_ij;
                achx_i -= ta_ij * dpx_ij;
                achy_i -= ta_ij * dpy_ij;
                achz_i -= ta_ij * dpz_ij;

                vndf vsu2_ij  = vndf::fabs(pp_i - vndf(epj[j].pres)) * rhi_ij * vndf(2.);
#ifdef __AVX512F__
                vndf vsui_ij  = _mm512_mask_mov_pd(vndf(0.),
                    _mm512_cmp_pd_mask(vsu2_ij, vndf(0.), _CMP_NEQ_OQ),
                    rsqrt(vsu2_ij));
#else
                vndf vsui_ij  = ((vsu2_ij != 0.) & rsqrt(vsu2_ij));
#endif
                vndf vsu_ij   = vsu2_ij * vsui_ij;
                vndf du_ij    = eg_i - vndf(epj[j].thrm);

                udot_i += ka_ij * (pc_i * w_ij
                                   - rhi_ij * (vndf(0.25) * av0_ij * w0_ij
                                               - (au_i + vndf(epj[j].alphu)) * vsu_ij * du_ij));

                difu_i += rhi_ij * du_ij * ka_ij * ri_ij;

                vndf dg_ij = vndf(epj[j].mass) * ri_ij * (et_i * dw_i + vndf(epj[j].eta) * dw_j);
                g1x_i += dg_ij * dpx_ij;
                g1y_i += dg_ij * dpy_ij;
                g1z_i += dg_ij * dpz_ij;
                    
            }

            achx_i *= vndf(0.5);
            achy_i *= vndf(0.5);
            achz_i *= vndf(0.5);
            udot_i *= vndf(0.5);
            difu_i *= vndf(2.);
            g1x_i  *= vndf(0.5);
            g1y_i  *= vndf(0.5);
            g1z_i  *= vndf(0.5);

            PS::F64 buf_ax[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_ay[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_az[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_du[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_vs[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_df[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_gx[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_gy[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_gz[nvector] __attribute__((aligned(vndf::nvector*8)));
            achx_i.store(buf_ax);
            achy_i.store(buf_ay);
            achz_i.store(buf_az);
            udot_i.store(buf_du);
            vsmx_i.store(buf_vs);
            difu_i.store(buf_df);
            g1x_i.store(buf_gx);
            g1y_i.store(buf_gy);
            g1z_i.store(buf_gz);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                hydro[i+ii].acch[0] = buf_ax[ii];
                hydro[i+ii].acch[1] = buf_ay[ii];
                hydro[i+ii].acch[2] = buf_az[ii];
                hydro[i+ii].udot    = buf_du[ii];
                hydro[i+ii].vsmx    = buf_vs[ii];
                hydro[i+ii].diffu   = buf_df[ii];
                hydro[i+ii].accg[0] = buf_gx[ii];
                hydro[i+ii].accg[1] = buf_gy[ii];
                hydro[i+ii].accg[2] = buf_gz[ii];
            }
            
        }

    }
};

#else

struct calcHydro {

    void operator () (const HydroEPI * epi,
                      const PS::S32 nip,
                      const HydroEPJ * epj,
                      const PS::S32 njp,
                      Hydro * hydro) {

        for(PS::S32 i = 0; i < nip; i++) {
            PS::S64 id_i = epi[i].id;
            PS::F64 px_i = epi[i].pos[0];
            PS::F64 py_i = epi[i].pos[1];
            PS::F64 pz_i = epi[i].pos[2];
            PS::F64 vx_i = epi[i].vel[0];
            PS::F64 vy_i = epi[i].vel[1];
            PS::F64 vz_i = epi[i].vel[2];
            PS::F64 hi_i = epi[i].hinv;
            PS::F64 pp_i = epi[i].pres;
            PS::F64 pc_i = epi[i].presc;
            PS::F64 bs_i = epi[i].bswt;
            PS::F64 dn_i = epi[i].dens;
            PS::F64 cs_i = epi[i].vsnd;
            PS::F64 al_i = epi[i].alph;
            PS::F64 au_i = epi[i].alphu;
            PS::F64 eg_i = epi[i].thrm;
            PS::F64 et_i = epi[i].eta;

            PS::F64 hi4_i  = hi_i * ND::calcVolumeInverse(hi_i);
            PS::F64 achx_i = 0.;
            PS::F64 achy_i = 0.;
            PS::F64 achz_i = 0.;
            PS::F64 udot_i = 0.;
            PS::F64 vsmx_i = 0.;
            PS::F64 difu_i = 0.;
            PS::F64 g1x_i  = 0.;
            PS::F64 g1y_i  = 0.;
            PS::F64 g1z_i  = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::F64 dpx_ij = px_i - epj[j].pos[0];
                PS::F64 dpy_ij = py_i - epj[j].pos[1];
                PS::F64 dpz_ij = pz_i - epj[j].pos[2];
                PS::F64 dvx_ij = vx_i - epj[j].vel[0];
                PS::F64 dvy_ij = vy_i - epj[j].vel[1];
                PS::F64 dvz_ij = vz_i - epj[j].vel[2];

                PS::F64 r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                PS::F64 ri_ij = ((id_i != epj[j].id) ? 1. / sqrt(r2_ij) : 0.);
                PS::F64 r1_ij = r2_ij * ri_ij;
                PS::F64 q_i   = r1_ij * hi_i;
                PS::F64 q_j   = r1_ij * epj[j].hinv;

                PS::F64 hi4_j = epj[j].hinv * ND::calcVolumeInverse(epj[j].hinv);
                PS::F64 dw_i  = hi4_i * SK::kernel1st(q_i);
                PS::F64 dw_j  = hi4_j * SK::kernel1st(q_j);
                PS::F64 ka_ij = (dw_i + dw_j) * epj[j].mass;

                PS::F64 rv_ij  = dpx_ij * dvx_ij + dpy_ij * dvy_ij + dpz_ij * dvz_ij;
                PS::F64 w_ij   = rv_ij * ri_ij;
                PS::F64 w0_ij  = ((w_ij < 0.) ? w_ij : 0.);
                PS::F64 vs_ij  = cs_i + epj[j].vsnd - 3. * w0_ij;
                PS::F64 rhi_ij = 1. / (dn_i + epj[j].dens);
                PS::F64 av0_ij = (bs_i + epj[j].bswt) * (al_i + epj[j].alph)
                    * vs_ij * w0_ij;
                                
                vsmx_i  = ((vsmx_i > vs_ij) ? vsmx_i : vs_ij);

                PS::F64 ta_ij = (pc_i + epj[j].presc - 0.5 * av0_ij * rhi_ij)
                    * ka_ij * ri_ij;
                achx_i -= ta_ij * dpx_ij;
                achy_i -= ta_ij * dpy_ij;
                achz_i -= ta_ij * dpz_ij;

                PS::F64 vsu2_ij = fabs(pp_i - epj[j].pres) * rhi_ij * 2.;
                PS::F64 vsui_ij = ((vsu2_ij != 0.) ? 1. / sqrt(vsu2_ij) : 0.);
                PS::F64 vsu_ij  = vsu2_ij * vsui_ij;
                PS::F64 du_ij   = eg_i - epj[j].thrm;

                udot_i += ka_ij * (pc_i * w_ij
                                   - rhi_ij * (0.25 * av0_ij * w0_ij
                                               - (au_i + epj[j].alphu) * vsu_ij * du_ij));

                difu_i += rhi_ij * du_ij * ka_ij * ri_ij;

                PS::F64 dg_ij = epj[j].mass * ri_ij * (et_i * dw_i + epj[j].eta * dw_j);
                g1x_i += dg_ij * dpx_ij;
                g1y_i += dg_ij * dpy_ij;
                g1z_i += dg_ij * dpz_ij;
                    
            }

            achx_i *= 0.5;
            achy_i *= 0.5;
            achz_i *= 0.5;
            udot_i *= 0.5;
            difu_i *= 2.;
            g1x_i  *= 0.5;
            g1y_i  *= 0.5;
            g1z_i  *= 0.5;

            hydro[i].acch[0] = achx_i;
            hydro[i].acch[1] = achy_i;
            hydro[i].acch[2] = achz_i;
            hydro[i].udot    = udot_i;
            hydro[i].vsmx    = vsmx_i;
            hydro[i].diffu   = difu_i;
            hydro[i].accg[0] = g1x_i;
            hydro[i].accg[1] = g1y_i;
            hydro[i].accg[2] = g1z_i;
            
        }

    }
};

#endif
