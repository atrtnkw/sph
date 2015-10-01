#pragma once

class Hydro {
public:
    PS::F64vec acch;
    PS::F64vec accg;
    PS::F64    udot;
    PS::F64    ydot;
    PS::F64    vsmx;
    void clear() {
        this->acch = 0.;
        this->accg = 0.;
        this->udot = 0.;
        this->ydot = 0.;
        this->vsmx = 0.;
    }
};

void SPH::copyFromForce(const Hydro & hydro) {
    this->acch  = hydro.acch;
    this->accg1 = hydro.accg;
    this->udot  = hydro.udot;
    this->ydot  = hydro.ydot;
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
    PS::F64    vary;
    PS::F64mat ctau;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;
    PS::F64    gamm;

    void copyFromFP(const SPH & sph){ 
        this->id   = sph.id;
        this->mass = sph.mass;
        this->pos  = sph.pos;
        this->vel  = sph.vel;
        this->ksr  = sph.ksr;
        this->hinv = 1. / sph.ksr;
        this->pres = sph.vary / sph.pres;
        this->vary = sph.vary;
        this->ctau = sph.ctau;
        this->bswt = 0.5 * sph.bswt;
        this->dens = sph.dens;
        this->vsnd = sph.vsnd;
        this->alph = sph.alph;
        this->gamm = RP::AdiabaticIndex;
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
    PS::F64    vary;
    PS::F64mat ctau;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;

    void copyFromFP(const SPH & sph){ 
        this->id   = sph.id;
        this->mass = sph.mass;
        this->pos  = sph.pos;
        this->vel  = sph.vel;
        this->ksr  = sph.ksr;
        this->hinv = 1. / sph.ksr;
        this->pres = sph.vary / sph.pres;
        this->vary = sph.vary;
        this->ctau = sph.ctau;
        this->bswt = 0.5 * sph.bswt;
        this->dens = sph.dens;
        this->vsnd = sph.vsnd;
        this->alph = sph.alph;
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
            PS::F64 buf_ys[nvector];
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
            PS::F64 buf_gm[nvector];
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
                buf_ys[ii] = epi[i+ii].vary;
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
                buf_gm[ii] = epi[i+ii].gamm - 1.;
            }
            v4df id_i, ms_i;
            v4df px_i, py_i, pz_i;
            v4df vx_i, vy_i, vz_i;
            v4df hi_i, ps_i, ys_i;
            v4df xx_i, yy_i, zz_i, xy_i, xz_i, yz_i;
            v4df bs_i, dn_i, cs_i, al_i;
            v4df gm_i;
            id_i.load(buf_id);
            ms_i.load(buf_ms);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            hi_i.load(buf_hi);
            ps_i.load(buf_ps);
            ys_i.load(buf_ys);
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
            gm_i.load(buf_gm);
            v4df hi3_i = ND::calcVolumeInverse(hi_i);;
            v4df achx_i(0.);
            v4df achy_i(0.);
            v4df achz_i(0.);
            v4df udot_i(0.);
            v4df ydot_i(0.);
            v4df vsmx_i(0.);
            for(PS::S32 j = 0; j < njp; j++) {
                v4df dpx_ij = px_i - v4df(epj[j].pos[0]);
                v4df dpy_ij = py_i - v4df(epj[j].pos[1]);
                v4df dpz_ij = pz_i - v4df(epj[j].pos[2]);
                v4df dvx_ij = vx_i - v4df(epj[j].vel[0]);
                v4df dvy_ij = vy_i - v4df(epj[j].vel[1]);
                v4df dvz_ij = vz_i - v4df(epj[j].vel[2]);

                v4df lax_i = dpx_ij * xx_i + dpy_ij * xy_i + dpz_ij * xz_i;
                v4df lay_i = dpx_ij * xy_i + dpy_ij * yy_i + dpz_ij * yz_i;
                v4df laz_i = dpx_ij * xz_i + dpy_ij * yz_i + dpz_ij * zz_i;
                v4df lax_j = dpx_ij * epj[j].ctau.xx + dpy_ij * epj[j].ctau.xy
                    + dpz_ij * epj[j].ctau.xz;
                v4df lay_j = dpx_ij * epj[j].ctau.xy + dpy_ij * epj[j].ctau.yy
                    + dpz_ij * epj[j].ctau.yz;
                v4df laz_j = dpx_ij * epj[j].ctau.xz + dpy_ij * epj[j].ctau.yz
                    + dpz_ij * epj[j].ctau.zz;

                v4df r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v4df ri_ij = ((id_i != v4df(epj[j].id)) & rsqrt(r2_ij));
                v4df r1_ij = r2_ij * ri_ij;
                v4df q_i   = r1_ij * hi_i;
                v4df q_j   = r1_ij * epj[j].hinv;

                v4df hi3_j = v4df(ND::calcVolumeInverse(epj[j].hinv));
                v4df w0_i  = v4df(-1.0) * hi3_i * SK::kernel0th(q_i);
                v4df w0_j  =              hi3_j * SK::kernel0th(q_j);

                lax_i *= w0_i;
                lay_i *= w0_i;
                laz_i *= w0_i;
                lax_j *= w0_j;
                lay_j *= w0_j;
                laz_j *= w0_j;                

                v4df lax_ij = v4df(0.5) * (lax_i - lax_j);
                v4df lay_ij = v4df(0.5) * (lay_i - lay_j);
                v4df laz_ij = v4df(0.5) * (laz_i - laz_j);

                v4df w_ij   = (dpx_ij * dvx_ij + dpy_ij * dvy_ij + dpz_ij * dvz_ij) * ri_ij;
                v4df w0_ij  = ((w_ij < v4df(0.)) & w_ij);
                v4df vs_ij  = cs_i + v4df(epj[j].vsnd) - v4df(3.) * w0_ij;
                v4df av_ij  = v4df(-0.5) * (bs_i + v4df(epj[j].bswt)) * (al_i + v4df(epj[j].alph))
                    * vs_ij * w0_ij * rcp(dn_i + v4df(epj[j].dens)) * ms_i * v4df(epj[j].mass);
                                
                vsmx_i  = v4df::max(vsmx_i, vs_ij);

                v4df ta_i = ps_i * v4df(epj[j].vary);
                lax_i  *= ta_i;
                lay_i  *= ta_i;
                laz_i  *= ta_i;
                lax_ij *= av_ij;
                lay_ij *= av_ij;
                laz_ij *= av_ij;

                achx_i -= (lax_i - v4df(epj[j].pres) * ys_i * lax_j + lax_ij);
                achy_i -= (lay_i - v4df(epj[j].pres) * ys_i * lay_j + lay_ij);
                achz_i -= (laz_i - v4df(epj[j].pres) * ys_i * laz_j + laz_ij);

                v4df vla_i  = dvx_ij * lax_i  + dvy_ij * lay_i  + dvz_ij * laz_i;
                v4df vla_ij = dvx_ij * lax_ij + dvy_ij * lay_ij + dvz_ij * laz_ij;
                udot_i += (vla_i + v4df(0.5) * vla_ij);

                ydot_i += vla_i;

            }

            v4df mi_i = rcp(ms_i);
            achx_i *= mi_i;
            achy_i *= mi_i;
            achz_i *= mi_i;
            udot_i *= mi_i;
            ydot_i *= gm_i;

            PS::F64 buf_ax[nvector];
            PS::F64 buf_ay[nvector];
            PS::F64 buf_az[nvector];
            PS::F64 buf_du[nvector];
            PS::F64 buf_dy[nvector];
            PS::F64 buf_vs[nvector];
            achx_i.store(buf_ax);
            achy_i.store(buf_ay);
            achz_i.store(buf_az);
            udot_i.store(buf_du);
            ydot_i.store(buf_dy);
            vsmx_i.store(buf_vs);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                hydro[i+ii].acch[0] = buf_ax[ii];
                hydro[i+ii].acch[1] = buf_ay[ii];
                hydro[i+ii].acch[2] = buf_az[ii];
                hydro[i+ii].udot    = buf_du[ii];
                hydro[i+ii].ydot    = buf_dy[ii];
                hydro[i+ii].vsmx    = buf_vs[ii];
            }
            
        }

    }
};
