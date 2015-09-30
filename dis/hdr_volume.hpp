#pragma once

class Volume {
public:
    PS::F64 vol;
    PS::F64 ksr;
    PS::F64 dens;
    PS::S64 np;
    void clear() {
        this->vol  = 0.0;
        this->ksr  = 0.0;
        this->dens = 0.0;
        this->np   = 0;
    }
};

void SPH::copyFromForce(const Volume & volume) {
    this->vol  = volume.vol;
    this->ksr  = volume.ksr;
    this->dens = volume.dens;
    this->np   = volume.np;
};

class VolumeEPI {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    ksr;
    PS::F64    pres;
    void copyFromFP(const SPH & sph){ 
        this->id   = sph.id;
        this->mass = sph.mass;
        this->pos  = sph.pos;
        this->ksr  = sph.ksr;
        this->pres = pow(sph.pres, 0.05);
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

class VolumeEPJ {
public:
    PS::S32    id;
    PS::F64vec pos;
    PS::F64    ksr;
    PS::F64    pres;
    void copyFromFP(const SPH & sph){ 
        this->id   = sph.id;
        this->pos  = sph.pos;
        this->ksr  = sph.ksr;
        this->pres = pow(sph.pres, 0.05);
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        this->pos = pos_new;
    }
};

struct calcVolume {

    void operator () (const VolumeEPI * epi,
                      const PS::S32 nip,
                      const VolumeEPJ * epj,
                      const PS::S32 njp,
                      Volume * volume) {

        v4df (*rcp)(v4df)   = v4df::rcp_4th;
        v4df (*rsqrt)(v4df) = v4df::rsqrt_4th;

        const PS::S32 nvector = v4df::getVectorLength();

        for(PS::S32 i = 0; i < nip; i += nvector) {
            const PS::S64 nii = std::min(nip - i, nvector);

            PS::F64 buf_id[nvector];
            PS::F64 buf_px[nvector];
            PS::F64 buf_py[nvector];
            PS::F64 buf_pz[nvector];
            PS::F64 buf_hs[nvector];
            PS::F64 buf_ps[nvector];
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_id[ii] = epi[i+ii].id;
                buf_px[ii] = epi[i+ii].pos[0];
                buf_py[ii] = epi[i+ii].pos[1];
                buf_pz[ii] = epi[i+ii].pos[2];
                buf_hs[ii] = epi[i+ii].ksr;
                buf_ps[ii] = epi[i+ii].pres;
            }
            v4df id_i;
            v4df px_i;
            v4df py_i;
            v4df pz_i;
            v4df hs_i;
            v4df ps_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            hs_i.load(buf_hs);
            ps_i.load(buf_ps);

            v4df hi_i  = rcp(hs_i);
            v4df hi3_i = ND::calcVolumeInverse(hi_i);
            v4df vl_i(0.);
            v4df nj_i(0.);

            for(PS::S32 j = 0; j < njp; j++) {
                v4df dx_ij = px_i - v4df(epj[j].pos[0]);
                v4df dy_ij = py_i - v4df(epj[j].pos[1]);
                v4df dz_ij = pz_i - v4df(epj[j].pos[2]);
                
                v4df r2_ij = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;                
                v4df r1_ij = v4df::sqrt(r2_ij);
                v4df q_i   = r1_ij * hi_i;                
                v4df w0_i  = hi3_i * SK::kernel0th(q_i);
                
                vl_i += v4df(epj[j].pres) * w0_i;
                nj_i += ((q_i < 1.) & v4df(1.));                    
            }

            vl_i = ps_i * rcp(vl_i);

            PS::F64 buf_vl[nvector], buf_nj[nvector];
            vl_i.store(buf_vl);
            nj_i.store(buf_nj);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                volume[i+ii].vol  = buf_vl[ii];
                volume[i+ii].ksr  = SK::ksrh * SK::eta * ND::calcPowerOfDimInverse(buf_vl[ii], 1.);
                volume[i+ii].dens = epi[i+ii].mass / buf_vl[ii];
                volume[i+ii].np   = buf_nj[ii];
            }

        }        
    }
};

template <class Tdinfo,
          class Tsph,
          class Tvolume>
void calcVolumeKernel(Tdinfo & dinfo,
                      Tsph & sph,
                      Tvolume & volume) {
    for(PS::S32 irepeat = 0; irepeat < 2; irepeat++) {
        volume.calcForceAllAndWriteBack(calcVolume(), sph, dinfo);        
        referEquationOfState(sph);
    }
}

