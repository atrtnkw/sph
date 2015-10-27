#pragma once

class Quantity {
public:
    PS::F64 stemp;
    PS::F64 scmps[NR::NumberOfNucleon];
    void clear() {
        this->stemp = 0.;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            this->scmps[k] = 0.;
        }
    }
};

void GeneralSPH::copyFromForce(const Quantity & quantity) {
    this->stemp = quantity.stemp;
    for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
        this->scmps[k] = quantity.scmps[k];
    }
}

class QuantityEPI {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64    hinv;
    PS::F64    ksr;

    void copyFromFP(const GeneralSPH & sph){ 
        this->id   = sph.id;
        this->pos  = sph.pos;
        this->hinv = 1. / sph.ksr;
        this->ksr  = sph.ksr;
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

class QuantityEPJ {
public:
    PS::S32    id;
    PS::F64    vol;
    PS::F64vec pos;
    PS::F64    hinv;
    PS::F64    hinv3;
    PS::F64    ksr;
    PS::F64    temp;
    PS::F64    cmps[NR::NumberOfNucleon];

    void copyFromFP(const GeneralSPH & sph){ 
        this->id    = sph.id;
        this->vol   = sph.mass / sph.dens;
        this->pos   = sph.pos;
        this->hinv  = 1. / sph.ksr;
        this->hinv3 = this->hinv * this->hinv * this->hinv;
        this->ksr   = sph.ksr;
        this->temp  = sph.temp;
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            this->cmps[k] = sph.cmps[k];
        }
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
            v4df hi_i(epi[i].hinv,   epi[i+1].hinv,   epi[i+2].hinv,   epi[i+3].hinv);
            v4df tt_i;
            v4df cm_i[NR::NumberOfNucleon];
            v4df hi3_i = hi_i * hi_i * hi_i;

            for(PS::S32 j = 0; j < njp; j++) {
                v4df dpx_ij = px_i - v4df(epj[j].pos[0]);
                v4df dpy_ij = py_i - v4df(epj[j].pos[1]);
                v4df dpz_ij = pz_i - v4df(epj[j].pos[2]);
                v4df r2_ij  = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v4df ri_ij  = ((v4df(epj[j].id) != id_i) & rsqrt(r2_ij));
                v4df r1_ij  = r2_ij * ri_ij;
                v4df q_i    = r1_ij * hi_i;
                v4df w_i    = hi3_i * SK::kernel0th(q_i);
                v4df q_j    = r1_ij * v4df(epj[j].hinv);
                v4df w_j    = v4df(epj[j].hinv3) * SK::kernel0th(q_j);
                v4df w_ij   = v4df(0.5) * (w_i + w_j);

                v4df vl_j(epj[j].vol);

                //tt_i += vl_j * v4df(epj[j].temp) * w_i;
                //tt_i += vl_j * v4df(epj[j].temp) * w_j;
                tt_i += vl_j * v4df(epj[j].temp) * w_ij;
                for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
                    //cm_i[k] += vl_j * v4df(epj[j].cmps[k]) * w_i;
                    //cm_i[k] += vl_j * v4df(epj[j].cmps[k]) * w_j;
                    cm_i[k] += vl_j * v4df(epj[j].cmps[k]) * w_ij;
                }
            }

            PS::F64 buf_temp[nvector];
            PS::F64 buf_cmps[NR::NumberOfNucleon][nvector];
            tt_i.store(buf_temp);
            for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
                cm_i[k].store(buf_cmps[k]);
            }

            PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;
            for(PS::S32 ii = 0; ii < nii; ii++) {
                quantity[i+ii].stemp = buf_temp[ii];
                for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
                    quantity[i+ii].scmps[k] = buf_cmps[k][ii];
                }
            }
        }

    }
};
