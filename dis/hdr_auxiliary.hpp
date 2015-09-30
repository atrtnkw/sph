#pragma once

class Auxiliary {
public:
    PS::F64    divv;
    PS::F64    rotv;
    PS::F64mat ctau;
    void clear() {
        this->divv = 0.0;
        this->rotv = 0.0;
        this->ctau = 0.0;
    }
};

void SPH::copyFromForce(const Auxiliary & auxiliary) {
    this->divv = auxiliary.divv;
    this->rotv = auxiliary.rotv;
    this->ctau = auxiliary.ctau;
};

class AuxiliaryEPI {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    pres;
    PS::F64    vol;
    void copyFromFP(const SPH & sph){ 
        this->id   = sph.id;
        this->pos  = sph.pos;
        this->vel  = sph.vel;
        this->ksr  = sph.ksr;
        this->pres = pow(sph.pres, 0.05);
        this->vol  = sph.vol;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        this->pos = pos_new;
    }
    PS::F64 getRSearch() const {
        return this->ksr;
    }
};

class AuxiliaryEPJ {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    pres;
    PS::F64    vol;
    void copyFromFP(const SPH & sph){ 
        this->id   = sph.id;
        this->pos  = sph.pos;
        this->vel  = sph.vel;
        this->ksr  = sph.ksr;
        this->pres = pow(sph.pres, 0.05);
        this->vol  = sph.vol;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        this->pos = pos_new;
    }
};

struct calcAuxiliary {

    void operator () (const AuxiliaryEPI * epi,
                      const PS::S32 nip,
                      const AuxiliaryEPJ * epj,
                      const PS::S32 njp,
                      Auxiliary * auxiliary) {

        v4df (*rcp)(v4df)   = v4df::rcp_4th;
        v4df (*rsqrt)(v4df) = v4df::rsqrt_4th;

        const PS::S32 nvector = v4df::getVectorLength();

        for(PS::S32 i = 0; i < nip; i += nvector) {
            const PS::S64 nii = std::min(nip - i, nvector);

            PS::F64 buf_id[nvector];
            PS::F64 buf_px[nvector];
            PS::F64 buf_py[nvector];
            PS::F64 buf_pz[nvector];
            PS::F64 buf_vx[nvector];
            PS::F64 buf_vy[nvector];
            PS::F64 buf_vz[nvector];
            PS::F64 buf_hs[nvector];
            PS::F64 buf_ps[nvector];
            PS::F64 buf_vl[nvector];
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_id[ii] = epi[i+ii].id;
                buf_px[ii] = epi[i+ii].pos[0];
                buf_py[ii] = epi[i+ii].pos[1];
                buf_pz[ii] = epi[i+ii].pos[2];
                buf_vx[ii] = epi[i+ii].vel[0];
                buf_vy[ii] = epi[i+ii].vel[1];
                buf_vz[ii] = epi[i+ii].vel[2];
                buf_hs[ii] = epi[i+ii].ksr;
                buf_ps[ii] = epi[i+ii].pres;
                buf_vl[ii] = epi[i+ii].vol;
            }
            v4df id_i;
            v4df px_i;
            v4df py_i;
            v4df pz_i;
            v4df vx_i;
            v4df vy_i;
            v4df vz_i;
            v4df hs_i;
            v4df ps_i;
            v4df vl_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            hs_i.load(buf_hs);
            ps_i.load(buf_ps);
            vl_i.load(buf_vl);

            v4df hi_i  = rcp(hs_i);
            v4df hi3_i = ND::calcVolumeInverse(hi_i);
            v4df hi4_i = hi_i * hi3_i;
            v4df divv_i(0.);
            v4df rotx_i(0.);
            v4df roty_i(0.);
            v4df rotz_i(0.);
            v4df ctxx_i(0.);
            v4df ctyy_i(0.);
            v4df ctzz_i(0.);
            v4df ctxy_i(0.);
            v4df ctxz_i(0.);
            v4df ctyz_i(0.);

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

                v4df w1_i = hi4_i * SK::kernel1st(q_i);
                v4df wc_i = ri_ij * w1_i * v4df(epj[j].pres);
                v4df wx_i = wc_i * dpx_ij;
                v4df wy_i = wc_i * dpy_ij;
                v4df wz_i = wc_i * dpz_ij;

                divv_i -= dvx_ij * wx_i;
                divv_i -= dvy_ij * wy_i;
                divv_i -= dvz_ij * wz_i;
                rotx_i += dvy_ij * wz_i - dvz_ij * wy_i;
                roty_i += dvz_ij * wx_i - dvx_ij * wz_i;
                rotz_i += dvx_ij * wy_i - dvy_ij * wx_i;

                v4df w0_i = hi3_i * SK::kernel0th(q_i);
                v4df vl_j = v4df(epj[j].vol) * w0_i;
                ctxx_i += vl_j * dpx_ij * dpx_ij;
                ctyy_i += vl_j * dpy_ij * dpy_ij;
                ctzz_i += vl_j * dpz_ij * dpz_ij;
                ctxy_i += vl_j * dpx_ij * dpy_ij;
                ctxz_i += vl_j * dpx_ij * dpz_ij;
                ctyz_i += vl_j * dpy_ij * dpz_ij;
            }            

            v4df rot2_i = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
            v4df rotv_i = rot2_i * ((rot2_i != 0.) & rsqrt(rot2_i));
            divv_i *= vl_i;
            rotv_i *= vl_i;

            PS::F64 buf_divv[nvector];
            PS::F64 buf_rotv[nvector];
            PS::F64 buf_ctxx[nvector];
            PS::F64 buf_ctyy[nvector];
            PS::F64 buf_ctzz[nvector];
            PS::F64 buf_ctxy[nvector];
            PS::F64 buf_ctxz[nvector];
            PS::F64 buf_ctyz[nvector];
            divv_i.store(buf_divv);
            rotv_i.store(buf_rotv);
            ctxx_i.store(buf_ctxx);
            ctyy_i.store(buf_ctyy);
            ctzz_i.store(buf_ctzz);
            ctxy_i.store(buf_ctxy);
            ctxz_i.store(buf_ctxz);
            ctyz_i.store(buf_ctyz);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                auxiliary[i+ii].divv = buf_divv[ii];
                auxiliary[i+ii].rotv = buf_rotv[ii];
                PS::F64mat otau = 0.;
                otau.xx = buf_ctxx[ii];
                otau.yy = buf_ctyy[ii];
                otau.zz = buf_ctzz[ii];
                otau.xy = buf_ctxy[ii];
                otau.xz = buf_ctxz[ii];
                otau.yz = buf_ctyz[ii];
                auxiliary[i+ii].ctau = ND::invertMatrix(otau);
            }

        }        
    }
};
