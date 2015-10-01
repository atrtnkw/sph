#pragma once

class Density {
public:
    PS::F64 dens;
    PS::F64 ksr;
    PS::S64 np;
    bool    itr;
    void clear() {
        dens = 0.;
        ksr  = 0.;
        np   = 0;
        itr  = false;
    }
};

void SPH::copyFromForce(const Density & density) {
    this->dens = density.dens;
    this->ksr  = density.ksr;
    this->np   = density.np;
};

class DensityEPI {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    ksr;
    PS::F64    rs;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
        ksr  = sph.ksr;
        rs   = sph.rs;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
    PS::F64 getRSearch() const {
        return this->rs;
    }
};

class DensityEPJ {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
};

struct calcDensityFirst {

    void operator () (const DensityEPI *epi,
                      const PS::S32 nip,
                      const DensityEPJ *epj,
                      const PS::S32 njp,
                      Density *density) {

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
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_id[ii] = epi[i+ii].id;
                buf_px[ii] = epi[i+ii].pos[0];
                buf_py[ii] = epi[i+ii].pos[1];
                buf_pz[ii] = epi[i+ii].pos[2];
                buf_hs[ii] = epi[i+ii].ksr;
            }
            v4df id_i;
            v4df px_i;
            v4df py_i;
            v4df pz_i;
            v4df hs_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            hs_i.load(buf_hs);

            for(PS::S32 repeat = 0; repeat < 3; repeat++) {
                v4df hi_i  = rcp(hs_i);
                v4df hi3_i = ND::calcVolumeInverse(hi_i);
                v4df hi4_i = hi_i * hi3_i;
                v4df rh_i(0.);
                v4df nj_i(0.);

                for(PS::S32 j = 0; j < njp; j++) {
                    v4df dx_ij = px_i - v4df(epj[j].pos[0]);
                    v4df dy_ij = py_i - v4df(epj[j].pos[1]);
                    v4df dz_ij = pz_i - v4df(epj[j].pos[2]);

                    v4df r2_ij = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;

                    v4df r1_ij = v4df::sqrt(r2_ij);
                    v4df q_i   = r1_ij * hi_i;

                    v4df kw0 = SK::kernel0th(q_i);
                    v4df rhj = v4df(epj[j].mass) * hi3_i * kw0;

                    rh_i += rhj;
                    nj_i += ((q_i < 1.) & v4df(1.));                    
                }

                PS::F64 buf_dens[nvector], buf_nj[nvector];
                rh_i.store(buf_dens);
                nj_i.store(buf_nj);
                for(PS::S32 ii = 0; ii < nii; ii++) {
                    buf_hs[ii] = SK::eta * SK::ksrh
                        * ND::calcPowerOfDimInverse(epi[i+ii].mass, buf_dens[ii]);
                    buf_hs[ii] = std::min(buf_hs[ii], RP::KernelSupportRadiusMaximum);
                    density[i+ii].dens = buf_dens[ii];
                    density[i+ii].np   = (PS::S64)buf_nj[ii];
                    density[i+ii].ksr  = buf_hs[ii];
                    density[i+ii].itr  = (buf_hs[ii] > epi[i+ii].rs) ? true : false;
                }
                hs_i.load(buf_hs);
            }
        }        
    }
};

template <class Tdinfo,
          class Tsph,
          class Tdensity>
void calcDensityKernel(Tdinfo & dinfo,
                       Tsph & sph,
                       Tdensity & density) {

    const PS::F64 expand  = 1.1;
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].rs = expand * sph[i].ksr;
    }

    PS::S32 cnt = 0;
    for(bool repeat = true; repeat == true;) {
        WT::start();
        bool repeat_loc = false;
        repeat = false;
        density.calcForceAll(calcDensityFirst(), sph, dinfo);
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].rs != 0.0) {
                if(density.getForce(i).itr == true) {
                    repeat_loc = true;
                    sph[i].rs *= expand;
                } else {
                    sph[i].rs = 0.0;
                    sph[i].copyFromForce(density.getForce(i));
                }
            }
        }
        WT::accumulateCalcDensity();
        repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
        cnt++;
    }
}
