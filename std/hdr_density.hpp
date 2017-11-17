#pragma once

class Density {
public:
    PS::F64 dens;
    PS::F64 divv;
    PS::F64 rotv;
    PS::F64 grdh;
    PS::F64 ksr;
    PS::S64 np;
    bool    itr;
    void clear() {
        dens = 0.;
        divv = 0.;
        rotv = 0.;
        grdh = 0.;
        ksr  = 0.;
        np   = 0;
        itr  = false;
    }
};

void SPH::copyFromForce(const Density & density) {
    this->dens = density.dens;
    this->divv = density.divv;
    this->rotv = density.rotv;
    this->grdh = density.grdh;
    this->ksr  = density.ksr;
    this->np   = density.np;
};

class DensityEPI {
public:
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    rs;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
        vel  = sph.vel;
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
    PS::F64vec vel;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
        vel  = sph.vel;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
};

#ifdef USE_INTRINSICS

struct calcDensity {

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

            PS::F64 buf_id[nvector] __attribute__((aligned(32)));
            PS::F64 buf_px[nvector] __attribute__((aligned(32)));
            PS::F64 buf_py[nvector] __attribute__((aligned(32)));
            PS::F64 buf_pz[nvector] __attribute__((aligned(32)));
            PS::F64 buf_vx[nvector] __attribute__((aligned(32)));
            PS::F64 buf_vy[nvector] __attribute__((aligned(32)));
            PS::F64 buf_vz[nvector] __attribute__((aligned(32)));
            PS::F64 buf_hs[nvector] __attribute__((aligned(32)));
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_id[ii] = epi[i+ii].id;
                buf_px[ii] = epi[i+ii].pos[0];
                buf_py[ii] = epi[i+ii].pos[1];
                buf_pz[ii] = epi[i+ii].pos[2];
                buf_vx[ii] = epi[i+ii].vel[0];
                buf_vy[ii] = epi[i+ii].vel[1];
                buf_vz[ii] = epi[i+ii].vel[2];
                buf_hs[ii] = epi[i+ii].ksr;
            }

            v4df id_i;
            v4df px_i;
            v4df py_i;
            v4df pz_i;
            v4df vx_i;
            v4df vy_i;
            v4df vz_i;
            v4df hs_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
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

                PS::F64 buf_dens[nvector] __attribute__((aligned(32)));
                PS::F64 buf_nj[nvector]   __attribute__((aligned(32)));
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
            
            v4df hi_i  = rcp(hs_i);
            v4df hi4_i = hi_i * ND::calcVolumeInverse(hi_i);
            v4df grdh_i(0.);
            v4df divv_i(0.);
            v4df rotx_i(0.);
            v4df roty_i(0.);
            v4df rotz_i(0.);
            for(PS::S32 j = 0; j < njp; j++) {
                v4df dpx_ij = px_i - v4df(epj[j].pos[0]);
                v4df dpy_ij = py_i - v4df(epj[j].pos[1]);
                v4df dpz_ij = pz_i - v4df(epj[j].pos[2]);
                v4df dvx_ij = vx_i - v4df(epj[j].vel[0]);
                v4df dvy_ij = vy_i - v4df(epj[j].vel[1]);
                v4df dvz_ij = vz_i - v4df(epj[j].vel[2]);

                v4df r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v4df ri_ij = rsqrt(r2_ij);
                ri_ij = ((id_i != v4df(epj[j].id)) & ri_ij);
                v4df r1_ij = r2_ij * ri_ij;
                v4df q_i = r1_ij * hi_i;

                v4df kw0 = SK::kernel0th(q_i);
                v4df kw1 = SK::kernel1st(q_i);

                v4df m_j(epj[j].mass);
                v4df ghj = v4df(RP::NumberOfDimension) * kw0;
                ghj += q_i * kw1;
                grdh_i -= ghj * hi4_i * m_j;

                v4df dw_ij  = m_j * hi4_i * kw1 * ri_ij;
                v4df dwx_ij = dw_ij * dpx_ij;
                v4df dwy_ij = dw_ij * dpy_ij;
                v4df dwz_ij = dw_ij * dpz_ij;

                divv_i -= dvx_ij * dwx_ij;
                divv_i -= dvy_ij * dwy_ij;
                divv_i -= dvz_ij * dwz_ij;

                rotx_i += dvy_ij * dwz_ij;
                roty_i += dvz_ij * dwx_ij;
                rotz_i += dvx_ij * dwy_ij;
                rotx_i -= dvz_ij * dwy_ij;
                roty_i -= dvx_ij * dwz_ij;
                rotz_i -= dvy_ij * dwx_ij;                 
            }

            v4df dens_i(density[i].dens, density[i+1].dens, density[i+2].dens, density[i+3].dens);
            v4df deni_i = rcp(dens_i);
            v4df omgi_i = rcp(v4df(1.) + hs_i * deni_i * grdh_i * rcp(RP::NumberOfDimension));
            v4df rot2_i = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
            v4df rotv_i = rot2_i * ((rot2_i != 0.) & rsqrt(rot2_i));
            rotv_i *= deni_i * omgi_i;
            divv_i *= deni_i * omgi_i;

            PS::F64 buf_divv[nvector] __attribute__((aligned(32)));
            PS::F64 buf_rotv[nvector] __attribute__((aligned(32)));
            PS::F64 buf_omgi[nvector] __attribute__((aligned(32)));
            divv_i.store(buf_divv);
            rotv_i.store(buf_rotv);
            omgi_i.store(buf_omgi);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                density[i+ii].divv = buf_divv[ii];
                density[i+ii].rotv = buf_rotv[ii];
                density[i+ii].grdh = buf_omgi[ii];
            }
        }        
    }
};

#else

struct calcDensity {

    void operator () (const DensityEPI *epi,
                      const PS::S32 nip,
                      const DensityEPJ *epj,
                      const PS::S32 njp,
                      Density *density) {

        for(PS::S32 i = 0; i < nip; i++) {
            PS::S64 id_i = epi[i].id;
            PS::F64 px_i = epi[i].pos[0];
            PS::F64 py_i = epi[i].pos[1];
            PS::F64 pz_i = epi[i].pos[2];
            PS::F64 vx_i = epi[i].vel[0];
            PS::F64 vy_i = epi[i].vel[1];
            PS::F64 vz_i = epi[i].vel[2];
            PS::F64 hs_i = epi[i].ksr;

            for(PS::S32 repeat = 0; repeat < 3; repeat++) {
                PS::F64 hi_i  = 1. / hs_i;
                PS::F64 hi3_i = ND::calcVolumeInverse(hi_i);
                PS::F64 hi4_i = hi_i * hi3_i;
                PS::F64 rh_i  = 0.;
                PS::S64 nj_i  = 0.;

                for(PS::S32 j = 0; j < njp; j++) {
                    PS::F64 dx_ij = px_i - epj[j].pos[0];
                    PS::F64 dy_ij = py_i - epj[j].pos[1];
                    PS::F64 dz_ij = pz_i - epj[j].pos[2];

                    PS::F64 r2_ij = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;

                    PS::F64 r1_ij = sqrt(r2_ij);
                    PS::F64 q_i   = r1_ij * hi_i;

                    PS::F64 kw0 = SK::kernel0th(q_i);
                    PS::F64 rhj = epj[j].mass * hi3_i * kw0;

                    rh_i += rhj;
                    nj_i += ((q_i < 1.) ? 1 : 0);                    
                }

                PS::F64 buf_hs = SK::eta * SK::ksrh
                    * ND::calcPowerOfDimInverse(epi[i].mass, rh_i);
                buf_hs = ((buf_hs < RP::KernelSupportRadiusMaximum)
                          ? buf_hs : RP::KernelSupportRadiusMaximum);
                density[i].dens = rh_i;
                density[i].np   = nj_i;
                density[i].ksr  = buf_hs;
                density[i].itr  = (buf_hs > epi[i].rs) ? true : false;

                hs_i = buf_hs;
            }
            
            PS::F64 hi_i   = 1. / hs_i;
            PS::F64 hi4_i  = hi_i * ND::calcVolumeInverse(hi_i);
            PS::F64 grdh_i = 0.;
            PS::F64 divv_i = 0.;
            PS::F64 rotx_i = 0.;
            PS::F64 roty_i = 0.;
            PS::F64 rotz_i = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::F64 dpx_ij = px_i - epj[j].pos[0];
                PS::F64 dpy_ij = py_i - epj[j].pos[1];
                PS::F64 dpz_ij = pz_i - epj[j].pos[2];
                PS::F64 dvx_ij = vx_i - epj[j].vel[0];
                PS::F64 dvy_ij = vy_i - epj[j].vel[1];
                PS::F64 dvz_ij = vz_i - epj[j].vel[2];

                PS::F64 r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                PS::F64 ri_ij = 1. / sqrt(r2_ij);
                ri_ij = ((id_i != epj[j].id) ? ri_ij : 0.);
                PS::F64 r1_ij = r2_ij * ri_ij;
                PS::F64 q_i = r1_ij * hi_i;

                PS::F64 kw0 = SK::kernel0th(q_i);
                PS::F64 kw1 = SK::kernel1st(q_i);

                PS::F64 m_j(epj[j].mass);
                PS::F64 ghj = PS::F64(RP::NumberOfDimension) * kw0;
                ghj += q_i * kw1;
                grdh_i -= ghj * hi4_i * m_j;

                PS::F64 dw_ij  = m_j * hi4_i * kw1 * ri_ij;
                PS::F64 dwx_ij = dw_ij * dpx_ij;
                PS::F64 dwy_ij = dw_ij * dpy_ij;
                PS::F64 dwz_ij = dw_ij * dpz_ij;

                divv_i -= dvx_ij * dwx_ij;
                divv_i -= dvy_ij * dwy_ij;
                divv_i -= dvz_ij * dwz_ij;

                rotx_i += dvy_ij * dwz_ij;
                roty_i += dvz_ij * dwx_ij;
                rotz_i += dvx_ij * dwy_ij;
                rotx_i -= dvz_ij * dwy_ij;
                roty_i -= dvx_ij * dwz_ij;
                rotz_i -= dvy_ij * dwx_ij;                 
            }

            PS::F64 dens_i = density[i].dens;
            PS::F64 deni_i = 1. / dens_i;
            PS::F64 omgi_i = 1. / (1. + hs_i * deni_i * grdh_i / RP::NumberOfDimension);
            PS::F64 rot2_i = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
            PS::F64 rotv_i = rot2_i * ((rot2_i != 0.) ? 1. / sqrt(rot2_i) : 0.);
            rotv_i *= deni_i * omgi_i;
            divv_i *= deni_i * omgi_i;

            density[i].divv = divv_i;
            density[i].rotv = rotv_i;
            density[i].grdh = omgi_i;
        }        
    }
};

#endif

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
        bool repeat_loc = false;
        repeat = false;
#ifdef DEBUG_REUSE_ON // a. tanikawa change here 17/11/17
        density.calcForceAll(calcDensity(), sph, dinfo, true, RP::getListMode());
#else
        density.calcForceAll(calcDensity(), sph, dinfo);
#endif
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
        repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
        cnt++;
    }
}
