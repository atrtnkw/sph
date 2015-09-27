#pragma once

class QuantityEPI {
public:
    PS::S32    id;
    PS::F64vec pos;

    void copyFromFP(const SPH & sph){ 
        id  = sph.id;
        pos = sph.pos;
    }

    void copyFromFP(const MassLess & msls){ 
        id  = msls.id;
        pos = msls.pos;
    }

    PS::F64vec getPos() const {
        return this->pos;
    }

    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }

};

class QuantityEPJ {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64    rhi;
    PS::F64vec pos;
    PS::F64    temp;
    PS::F64    abar;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    shck;

    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        rhi  = 1. / sph.dens;
        pos  = sph.pos;
        temp = sph.temp;
        abar = sph.abar;
        hinv = 1. / sph.ksr;
        ksr  = sph.ksr;
        shck = sph.bswt * sph.ksr * sph.divv / sph.vsnd;
    }

    void copyFromFP(const MassLess & msls){ 
        id   = msls.id;
        mass = msls.mass;
        rhi  = 0.;
        pos  = msls.pos;
        temp = 0.;
        abar = 0.;
        hinv = 0.;
        ksr  = 0.;
        shck = 0.;
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

struct calcQuantity {

    void operator () (const QuantityEPI *epi,
                      const PS::S32 nip,
                      const QuantityEPJ *epj,
                      const PS::S32 njp,
                      Quantity *quantity) {
        v4df (*rcp)(v4df)   = v4df::rcp_4th;
        v4df (*rsqrt)(v4df) = v4df::rsqrt_4th;

        PS::S32 nvector = v4df::getVectorLength();

        for(PS::S32 i = 0; i < nip; i += nvector) {
            v4df id_i(epi[i].id,     epi[i+1].id,      epi[i+2].id,    epi[i+3].id);
            v4df px_i(epi[i].pos[0], epi[i+1].pos[0], epi[i+2].pos[0], epi[i+3].pos[0]);
            v4df py_i(epi[i].pos[1], epi[i+1].pos[1], epi[i+2].pos[1], epi[i+3].pos[1]);
            v4df pz_i(epi[i].pos[2], epi[i+1].pos[2], epi[i+2].pos[2], epi[i+3].pos[2]);
            v4df dens_i(0.d);
            v4df temp_i(0.d);
            v4df abar_i(0.d);
            v4df shck_i(0.d);

            for(PS::S32 j = 0; j < njp; j++) {
                v4df dpx_ij = px_i - v4df(epj[j].pos[0]);
                v4df dpy_ij = py_i - v4df(epj[j].pos[1]);
                v4df dpz_ij = pz_i - v4df(epj[j].pos[2]);
                v4df m_j    = v4df(epj[j].mass);
                v4df hi_j   = v4df(epj[j].hinv);
                v4df hi3_j  = SPH::calcVolumeInverse(hi_j);

                v4df r2_ij  = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v4df ri_ij  = ((v4df(epj[j].id) != id_i) & rsqrt(r2_ij));
                v4df r1_ij  = r2_ij * ri_ij;
                v4df q_j    = r1_ij * hi_j;
                v4df w_j    = hi3_j * KernelSph::kernel0th(q_j);

                dens_i += m_j * w_j;
                temp_i += m_j * v4df(epj[j].rhi) * v4df(epj[j].temp) * w_j;
                abar_i += m_j * v4df(epj[j].rhi) * v4df(epj[j].abar) * w_j;
                shck_i += m_j * v4df(epj[j].rhi) * v4df(epj[j].shck) * w_j;

            }

            PS::F64 buf_dens[nvector];
            PS::F64 buf_temp[nvector];
            PS::F64 buf_abar[nvector];
            PS::F64 buf_shck[nvector];
            dens_i.store(buf_dens);
            temp_i.store(buf_temp);
            abar_i.store(buf_abar);
            shck_i.store(buf_shck);

            PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;
            for(PS::S32 ii = 0; ii < nii; ii++) {
                quantity[i+ii].dens = buf_dens[ii];
                quantity[i+ii].temp = buf_temp[ii];
                quantity[i+ii].abar = buf_abar[ii];
                quantity[i+ii].shck = buf_shck[ii];
            }

        }
    }
};
