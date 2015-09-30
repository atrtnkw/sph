#pragma once

class Hydro {
public:
    PS::F64vec acch;
    PS::F64vec accg;
    PS::F64    udot;
    PS::F64    vsmx;
    void clear() {
        this->acch = 0.;
        this->accg = 0.;
        this->udot = 0.;
        this->vsmx = 0.;
    }
};

void SPH::copyFromForce(const Hydro & hydro) {
    this->acch  = hydro.acch;
    this->accg1 = hydro.accg;
    this->udot  = hydro.udot;
    this->vsmx  = hydro.vsmx;
    this->acc   = this->acch + this->accg1 + this->accg2;
}

class HydroEPI {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    pres;
    PS::F64    presk;
    PS::F64    vol;
    PS::F64mat ctau;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;

    void copyFromFP(const SPH & sph){ 
        this->id    = sph.id;
        this->mass  = sph.mass;
        this->pos   = sph.pos;
        this->vel   = sph.vel;
        this->ksr   = sph.ksr;
        this->hinv  = 1. / sph.ksr;
        this->pres  = sph.pres;
        this->presk = pow(sph.pres, 0.05);
        this->vol   = sph.vol;
        this->ctau  = sph.ctau;
        this->bswt  = 0.5 * sph.bswt;
        this->dens  = sph.dens;
        this->vsnd  = sph.vsnd;
        this->alph  = sph.alph;
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
    PS::F64    presk;
    PS::F64    vol;
    PS::F64mat ctau;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;

    void copyFromFP(const SPH & sph){ 
        this->id    = sph.id;
        this->mass  = sph.mass;
        this->pos   = sph.pos;
        this->vel   = sph.vel;
        this->ksr   = sph.ksr;
        this->hinv  = 1. / sph.ksr;
        this->pres  = sph.pres;
        this->presk = pow(sph.pres, 0.05);
        this->vol   = sph.vol;
        this->ctau  = sph.ctau;
        this->bswt  = 0.5 * sph.bswt;
        this->dens  = sph.dens;
        this->vsnd  = sph.vsnd;
        this->alph  = sph.alph;
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

struct calcHydro {

    void operator () (const HydroEPI * epi,
                      const PS::S32 nip,
                      const HydroEPJ * epj,
                      const PS::S32 njp,
                      Hydro * hydro) {
        v4df (*rcp)(v4df)   = v4df::rcp_4th;
        v4df (*rsqrt)(v4df) = v4df::rsqrt_4th;

        PS::S32 nvector = v4df::getVectorLength();

        for(PS::S32 i = 0; i < nip; i += nvector) {
            PS::F64 buf_id[nvector];
            PS::F64 buf_ms[nvector];
            PS::F64 buf_px[nvector];
            PS::F64 buf_py[nvector];
            PS::F64 buf_pz[nvector];
            PS::F64 buf_vx[nvector];
            PS::F64 buf_vy[nvector];
            PS::F64 buf_vz[nvector];
            PS::F64 buf_hi[nvector];
            PS::F64 buf_ps[nvector];
            PS::F64 buf_pk[nvector];
            PS::F64 buf_vl[nvector];
            PS::F64 buf_xx[nvector];
            PS::F64 buf_yy[nvector];
            PS::F64 buf_zz[nvector];
            PS::F64 buf_xy[nvector];
            PS::F64 buf_xz[nvector];
            PS::F64 buf_yz[nvector];
            PS::F64 buf_bs[nvector];
            PS::F64 buf_dn[nvector];
            PS::F64 buf_cs[nvector];
            PS::F64 buf_al[nvector];
            PS::S32 nii = std::min(nip - i, nvector);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_id[ii] = epi[i+ii].id;
                buf_ms[ii] = epi[i+ii].mass;
                buf_px[ii] = epi[i+ii].pos[0];
                buf_py[ii] = epi[i+ii].pos[1];
                buf_pz[ii] = epi[i+ii].pos[2];
                buf_vx[ii] = epi[i+ii].vel[0];
                buf_vy[ii] = epi[i+ii].vel[1];
                buf_vz[ii] = epi[i+ii].vel[2];
                buf_hi[ii] = epi[i+ii].hinv;
                buf_ps[ii] = epi[i+ii].pres;
                buf_pk[ii] = epi[i+ii].presk;
                buf_vl[ii] = epi[i+ii].vol;
                buf_xx[ii] = epi[i+ii].ctau.xx;
                buf_yy[ii] = epi[i+ii].ctau.yy;
                buf_zz[ii] = epi[i+ii].ctau.zz;
                buf_xy[ii] = epi[i+ii].ctau.xy;
                buf_xz[ii] = epi[i+ii].ctau.xz;
                buf_yz[ii] = epi[i+ii].ctau.yz;
                buf_bs[ii] = epi[i+ii].bswt;
                buf_dn[ii] = epi[i+ii].dens;
                buf_cs[ii] = epi[i+ii].vsnd;
                buf_al[ii] = epi[i+ii].alph;
            }
            v4df id_i;
            v4df px_i, py_i, pz_i;
            v4df vx_i, vy_i, vz_i;
            v4df hi_i, vl_i;
            v4df ps_i, pk_i;            
            v4df xx_i, yy_i, zz_i, xy_i, xz_i, yz_i;
            v4df bs_i, dn_i, cs_i, al_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            hi_i.load(buf_hi);
            vl_i.load(buf_vl);
            ps_i.load(buf_ps);
            pk_i.load(buf_pk);
            xx_i.load(buf_xx);
            yy_i.load(buf_yy);
            zz_i.load(buf_zz);
            xy_i.load(buf_xy);
            xz_i.load(buf_xz);
            yz_i.load(buf_yz);
            bs_i.load(buf_bs);
            dn_i.load(buf_dn);
            cs_i.load(buf_cs);
            al_i.load(buf_al);
            v4df hi3_i = ND::calcVolumeInverse(hi_i);;
            v4df vl2_i = vl_i * vl_i;
            v4df achx_i(0.);
            v4df achy_i(0.);
            v4df achz_i(0.);
            v4df udot_i(0.);
            v4df vsmx_i(0.);
            for(PS::S32 j = 0; j < njp; j++) {
                v4df dpx_ij = px_i - v4df(epj[j].pos[0]);
                v4df dpy_ij = py_i - v4df(epj[j].pos[1]);
                v4df dpz_ij = pz_i - v4df(epj[j].pos[2]);
                v4df dvx_ij = vx_i - v4df(epj[j].vel[0]);
                v4df dvy_ij = vy_i - v4df(epj[j].vel[1]);
                v4df dvz_ij = vz_i - v4df(epj[j].vel[2]);

                v4df r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v4df ri_ij = ((id_i != v4df(epj[j].id)) & rsqrt(r2_ij));
                v4df r1_ij = r2_ij * ri_ij;
                v4df q_i   = r1_ij * hi_i;
                v4df q_j   = r1_ij * epj[j].hinv;

                v4df lax_i = dpx_ij * xx_i + dpy_ij * xy_i + dpz_ij * xz_i;
                v4df lay_i = dpx_ij * xy_i + dpy_ij * yy_i + dpz_ij * yz_i;
                v4df laz_i = dpx_ij * xz_i + dpy_ij * yz_i + dpz_ij * zz_i;
                v4df lax_j = dpx_ij * epj[j].ctau.xx + dpy_ij * epj[j].ctau.xy
                    + dpz_ij * epj[j].ctau.xz;
                v4df lay_j = dpx_ij * epj[j].ctau.xy + dpy_ij * epj[j].ctau.yy
                    + dpz_ij * epj[j].ctau.yz;
                v4df laz_j = dpx_ij * epj[j].ctau.xz + dpy_ij * epj[j].ctau.yz
                    + dpz_ij * epj[j].ctau.zz;
                
                v4df hi3_j = v4df(ND::calcVolumeInverse(epj[j].hinv));
                v4df w0_i  = v4df(-1.0) * hi3_i * SK::kernel0th(q_i);
                v4df w0_j  = hi3_j * SK::kernel0th(q_j);

                lax_i *= w0_i;
                lay_i *= w0_i;
                laz_i *= w0_i;
                lax_j *= w0_j;
                lay_j *= w0_j;
                laz_j *= w0_j;

                v4df lax_ij = v4df(0.5) * (lax_i - lax_j);
                v4df lay_ij = v4df(0.5) * (lay_i - lay_j);
                v4df laz_ij = v4df(0.5) * (laz_i - laz_j);

                v4df cax_i  = ps_i * vl2_i * ;

            }
        }
        /*
        for(PS::S32 i = 0; i < nip; i += nvector) {
            for(PS::S32 j = 0; j < njp; j++) {
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

                v4df de_ij = prhi2_i + v4df(0.5d) * vis_ij;
                de_ij *= mrdw_ij;
                ene_i += xv_ij * de_ij;                

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
        */

    }
};
