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
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    presc;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;

    void copyFromFP(const SPH & sph){ 
        this->id    = sph.id;
        this->pos   = sph.pos;
        this->vel   = sph.vel;
        this->ksr   = sph.ksr;
        this->hinv  = 1. / sph.ksr;
        this->presc = sph.pres * sph.grdh / (sph.dens * sph.dens);
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
    PS::F64    presc;
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
        this->presc = sph.pres * sph.grdh / (sph.dens * sph.dens);
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
            PS::F64 buf_px[nvector];
            PS::F64 buf_py[nvector];
            PS::F64 buf_pz[nvector];
            PS::F64 buf_vx[nvector];
            PS::F64 buf_vy[nvector];
            PS::F64 buf_vz[nvector];
            PS::F64 buf_hi[nvector];
            PS::F64 buf_pc[nvector];
            PS::F64 buf_bs[nvector];
            PS::F64 buf_dn[nvector];
            PS::F64 buf_cs[nvector];
            PS::F64 buf_al[nvector];
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
                buf_pc[ii] = epi[i+ii].presc;
                buf_bs[ii] = epi[i+ii].bswt;
                buf_dn[ii] = epi[i+ii].dens;
                buf_cs[ii] = epi[i+ii].vsnd;
                buf_al[ii] = epi[i+ii].alph;
            }
            v4df id_i;
            v4df px_i, py_i, pz_i;
            v4df vx_i, vy_i, vz_i;
            v4df hi_i, pc_i;
            v4df bs_i, dn_i, cs_i, al_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            hi_i.load(buf_hi);
            pc_i.load(buf_pc);
            bs_i.load(buf_bs);
            dn_i.load(buf_dn);
            cs_i.load(buf_cs);
            al_i.load(buf_al);
            v4df hi4_i = hi_i * ND::calcVolumeInverse(hi_i);
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

                v4df hi4_j = v4df(epj[j].hinv) * ND::calcVolumeInverse(epj[j].hinv);
                v4df ka_ij = (hi4_i * SK::kernel1st(q_i) + hi4_j * SK::kernel1st(q_j))
                    * ri_ij * v4df(epj[j].mass);

                v4df rv_ij = dpx_ij * dvx_ij + dpy_ij * dvy_ij + dpz_ij * dvz_ij;
                v4df w_ij  = rv_ij * ri_ij;
                v4df w0_ij = ((w_ij < v4df(0.)) & w_ij);
                v4df vs_ij = cs_i + v4df(epj[j].vsnd) - v4df(3.) * w0_ij;
                v4df av_ij = v4df(-0.5) * (bs_i + v4df(epj[j].bswt)) * (al_i + v4df(epj[j].alph))
                    * vs_ij * w0_ij * rcp(dn_i + v4df(epj[j].dens));
                                
                vsmx_i  = v4df::max(vsmx_i, vs_ij);

                v4df ta_ij = (pc_i + v4df(epj[j].presc) + av_ij) * ka_ij;
                achx_i -= ta_ij * dpx_ij;
                achy_i -= ta_ij * dpy_ij;
                achz_i -= ta_ij * dpz_ij;
                udot_i += (pc_i + v4df(0.5) * av_ij) * ka_ij * rv_ij;
                    
            }

            achx_i *= v4df(0.5);
            achy_i *= v4df(0.5);
            achz_i *= v4df(0.5);
            udot_i *= v4df(0.5);

            PS::F64 buf_ax[nvector];
            PS::F64 buf_ay[nvector];
            PS::F64 buf_az[nvector];
            PS::F64 buf_du[nvector];
            PS::F64 buf_vs[nvector];
            achx_i.store(buf_ax);
            achy_i.store(buf_ay);
            achz_i.store(buf_az);
            udot_i.store(buf_du);
            vsmx_i.store(buf_vs);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                hydro[i+ii].acch[0] = buf_ax[ii];
                hydro[i+ii].acch[1] = buf_ay[ii];
                hydro[i+ii].acch[2] = buf_az[ii];
                hydro[i+ii].udot    = buf_du[ii];
                hydro[i+ii].vsmx    = buf_vs[ii];
            }
            
        }

    }
};
