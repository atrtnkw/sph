#pragma once

class Auxiliary {
public:
    PS::F64    divv;
    PS::F64    rotv;
    PS::F64    dens;
    PS::F64    grdh;
    PS::S64    nj2;
    void clear() {
        this->divv = 0.0;
        this->rotv = 0.0;
        this->dens = 0.0;
        this->grdh = 0.0;
        this->nj2  = 0;
    }
};

void SPH::copyFromForce(const Auxiliary & auxiliary) {
    this->divv = auxiliary.divv;
    this->rotv = auxiliary.rotv;
    this->dens = auxiliary.dens;
    //this->grdh = auxiliary.grdh;
    this->grdh = 1.;
    this->np2  = auxiliary.nj2;
};

class AuxiliaryEPI {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    pres;
    void copyFromFP(const SPH & sph){ 
        this->id   = sph.id;
        this->pos  = sph.pos;
        this->vel  = sph.vel;
        this->ksr  = sph.ksr;
        this->pres = sph.presk;
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
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    pres;
    PS::F64    vary;
    void copyFromFP(const SPH & sph){ 
        this->id   = sph.id;
        this->mass = sph.mass;
        this->pos  = sph.pos;
        this->vel  = sph.vel;
        this->pres = sph.presk;
        this->vary = sph.vary;
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
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            hs_i.load(buf_hs);
            ps_i.load(buf_ps);

            v4df hi_i  = rcp(hs_i);
            v4df hi3_i = ND::calcVolumeInverse(hi_i);
            v4df hi4_i = hi_i * hi3_i;
            v4df divv_i(0.);
            v4df rotx_i(0.);
            v4df roty_i(0.);
            v4df rotz_i(0.);
            v4df dens_i(0.);
            v4df grdh_i(0.);
            v4df nj2_i(0.);

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

                v4df w1_i = SK::kernel1st(q_i);
                v4df wc_i = ri_ij * hi4_i * w1_i * v4df(epj[j].vary);
                v4df wx_i = wc_i * dpx_ij;
                v4df wy_i = wc_i * dpy_ij;
                v4df wz_i = wc_i * dpz_ij;

                divv_i -= dvx_ij * wx_i;
                divv_i -= dvy_ij * wy_i;
                divv_i -= dvz_ij * wz_i;
                rotx_i += dvy_ij * wz_i - dvz_ij * wy_i;
                roty_i += dvz_ij * wx_i - dvx_ij * wz_i;
                rotz_i += dvx_ij * wy_i - dvy_ij * wx_i;

                v4df w0_i = SK::kernel0th(q_i);
                dens_i += v4df(epj[j].mass) * hi3_i * w0_i;
                grdh_i -= v4df(epj[j].vary) * hi4_i * (w0_i * RP::NumberOfDimension + q_i * w1_i);

                nj2_i  += ((q_i < 1.) & v4df(1.));
            }            

            v4df omgi_i = rcp(v4df(1.) + hs_i * grdh_i * rcp(ps_i * RP::NumberOfDimension));
            v4df rot2_i = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
            v4df rotv_i = rot2_i * ((rot2_i != 0.) & rsqrt(rot2_i));
            //v4df pi_i = v4df(-1.) * rcp(ps_i);
            v4df pi_i   = v4df(-1.) * rcp(ps_i) * omgi_i;
            divv_i *= pi_i;
            rotv_i *= pi_i;

            PS::F64 buf_divv[nvector];
            PS::F64 buf_rotv[nvector];
            PS::F64 buf_dens[nvector];
            PS::F64 buf_omgi[nvector];
            PS::F64 buf_nj2[nvector];
            divv_i.store(buf_divv);
            rotv_i.store(buf_rotv);
            dens_i.store(buf_dens);
            omgi_i.store(buf_omgi);
            nj2_i.store(buf_nj2);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                auxiliary[i+ii].divv = buf_divv[ii];
                auxiliary[i+ii].rotv = buf_rotv[ii];
                auxiliary[i+ii].dens = buf_dens[ii];
                auxiliary[i+ii].grdh = buf_omgi[ii];
                auxiliary[i+ii].nj2  = buf_nj2[ii];
            }

        }        
    }
};
