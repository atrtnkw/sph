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

//#define USE_SIMDI4J2
#ifdef USE_SIMDI4J2

struct calcDensity {

    static const PS::S64 njarraymax = 65536;

    void operator () (const DensityEPI *epi,
                      const PS::S32 nip,
                      const DensityEPJ *epj,
                      const PS::S32 njp,
                      Density *density) {

        static __thread PS::S32 idj[njarraymax];
        static __thread PS::F64 mj[njarraymax];
        static __thread PS::F64 pxj[njarraymax];
        static __thread PS::F64 pyj[njarraymax];
        static __thread PS::F64 pzj[njarraymax];
        static __thread PS::F64 vxj[njarraymax];
        static __thread PS::F64 vyj[njarraymax];
        static __thread PS::F64 vzj[njarraymax];

        v8df (*rcp)(v8df)   = v8df::rcp_4th;
        v8df (*rsqrt)(v8df) = v8df::rsqrt_4th;

        const PS::S32 nivec   = 4;
        const PS::S32 njvec   = 2;
        const PS::S32 njarray = (njp / njvec + 1) * njvec;
        assert(njarray <= njarraymax);
        for(PS::S32 j = 0; j < njarray; j += njvec) {
            for(PS::S32 jj = 0; jj < njvec; jj++) {
                idj[j+jj] = epj[j+jj].id;
                mj[j+jj]  = epj[j+jj].mass;
                pxj[j+jj] = epj[j+jj].pos[0];
                pyj[j+jj] = epj[j+jj].pos[1];
                pzj[j+jj] = epj[j+jj].pos[2];
                vxj[j+jj] = epj[j+jj].vel[0];
                vyj[j+jj] = epj[j+jj].vel[1];
                vzj[j+jj] = epj[j+jj].vel[2];
            }
        }
        if(njarray > njp) {
            idj[njarray-1]      = -1;
            mj[njarray-1]       = 0.;
            pxj[njarray-1]      = 1e+30;
        }

        for(PS::S32 i = 0; i < nip; i += nivec) {
            const PS::S32 nii = std::min(nip - i, nivec);
            PS::F64 buf_id[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_px[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_py[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_pz[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_vx[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_vy[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_vz[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_hs[nivec] __attribute__((aligned(nivec*8)));
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
                v8df hi_i  = rcp(v8df(hs_i));
                v8df hi3_i = ND::calcVolumeInverse(hi_i);
                v8df hi4_i = hi_i * hi3_i;
                v8df rh_i(0.);
                v8df nj_i(0.);
                for(PS::S32 j = 0; j < njarray; j += njvec) {
                    v8df dx_ij = v8df(px_i) - v8df(v4df(pxj[j]), v4df(pxj[j+1]));
                    v8df dy_ij = v8df(py_i) - v8df(v4df(pyj[j]), v4df(pyj[j+1]));
                    v8df dz_ij = v8df(pz_i) - v8df(v4df(pzj[j]), v4df(pzj[j+1]));
                    v8df r2_ij = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;
                    v8df r1_ij = v8df::sqrt(r2_ij);
                    v8df q_i   = r1_ij * hi_i;
                    v8df kw0   = SK::kernel0th(q_i);
                    v8df rhj   = v8df(v4df(mj[j]), v4df(mj[j+1])) * hi3_i * kw0;
                    rh_i += rhj;
                    nj_i += _mm512_mask_mov_pd(v8df(0.),
                                               _mm512_cmp_pd_mask(q_i, v8df(1.), _CMP_LT_OS),
                                               v8df(1.));
                }
                PS::F64 buf_dens[nivec] __attribute__((aligned(nivec*8)));
                PS::F64 buf_nj[nivec]   __attribute__((aligned(nivec*8)));
                v8df::reduce(rh_i).store(buf_dens);
                v8df::reduce(nj_i).store(buf_nj);
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

            v8df hi_i  = rcp(v8df(hs_i));
            v8df hi4_i = hi_i * ND::calcVolumeInverse(hi_i);
            v8df grdh_i(0.);
            v8df divv_i(0.);
            v8df rotx_i(0.);
            v8df roty_i(0.);
            v8df rotz_i(0.);
            for(PS::S32 j = 0; j < njarray; j += njvec) {
                v8df dpx_ij = v8df(px_i) - v8df(v4df(pxj[j]), v4df(pxj[j+1]));
                v8df dpy_ij = v8df(py_i) - v8df(v4df(pyj[j]), v4df(pyj[j+1]));
                v8df dpz_ij = v8df(pz_i) - v8df(v4df(pzj[j]), v4df(pzj[j+1]));
                v8df dvx_ij = v8df(vx_i) - v8df(v4df(vxj[j]), v4df(vxj[j+1]));
                v8df dvy_ij = v8df(vy_i) - v8df(v4df(vyj[j]), v4df(vyj[j+1]));
                v8df dvz_ij = v8df(vz_i) - v8df(v4df(vzj[j]), v4df(vzj[j+1]));

                v8df r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v8df ri_ij = rsqrt(r2_ij);
                ri_ij = _mm512_mask_mov_pd(v8df(0.),
                                           _mm512_cmp_pd_mask(v8df(id_i),
                                                              v8df(v4df(idj[j]), v4df(idj[j+1])),
                                                              _CMP_NEQ_UQ),
                                           ri_ij);
                v8df r1_ij = r2_ij * ri_ij;
                v8df q_i = r1_ij * hi_i;

                v8df kw0 = SK::kernel0th(q_i);
                v8df kw1 = SK::kernel1st(q_i);

                v8df m_j(v4df(mj[j]), v4df(mj[j+1]));
                v8df ghj = v8df(RP::NumberOfDimension) * kw0;
                ghj += q_i * kw1;
                grdh_i -= ghj * hi4_i * m_j;

                v8df dw_ij  = m_j * hi4_i * kw1 * ri_ij;
                v8df dwx_ij = dw_ij * dpx_ij;
                v8df dwy_ij = dw_ij * dpy_ij;
                v8df dwz_ij = dw_ij * dpz_ij;

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

            v4df grdh = v8df::reduce(grdh_i);
            v4df divv = v8df::reduce(divv_i);
            v4df rotx = v8df::reduce(rotx_i);
            v4df roty = v8df::reduce(roty_i);
            v4df rotz = v8df::reduce(rotz_i);
            PS::F64 buf_dens[nivec] __attribute__((aligned(nivec*8)));
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_dens[ii] = density[i+ii].dens;
            }
            v4df dens;
            dens.load(buf_dens);
            v4df deni = v4df::rcp_4th(dens);
            v4df omgi = v4df::rcp_4th(v4df(1.) + hs_i * deni * grdh
                                      * v4df::rcp_4th(RP::NumberOfDimension));
            v4df rot2 = rotx * rotx + roty * roty + rotz * rotz;
            v4df rotv = rot2 * _mm256_mask_mov_pd(v4df(0.),
                                                  _mm256_cmp_pd_mask(rot2, v4df(0.), _CMP_NEQ_UQ),
                                                  v4df::rsqrt_4th(rot2));
            rotv *= deni * omgi;
            divv *= deni * omgi;
            PS::F64 buf_divv[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_rotv[nivec] __attribute__((aligned(nivec*8)));
            PS::F64 buf_omgi[nivec] __attribute__((aligned(nivec*8)));
            divv.store(buf_divv);
            rotv.store(buf_rotv);
            omgi.store(buf_omgi);
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

        vndf (*rcp)(vndf)   = vndf::rcp_4th;
        vndf (*rsqrt)(vndf) = vndf::rsqrt_4th;

        const PS::S32 nvector = vndf::nvector;
        
        for(PS::S32 i = 0; i < nip; i += nvector) {
            const PS::S64 nii = std::min(nip - i, nvector);

            PS::F64 buf_id[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_px[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_py[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_pz[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_vx[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_vy[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_vz[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_hs[nvector] __attribute__((aligned(vndf::nvector*8)));
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

            vndf id_i;
            vndf px_i;
            vndf py_i;
            vndf pz_i;
            vndf vx_i;
            vndf vy_i;
            vndf vz_i;
            vndf hs_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            hs_i.load(buf_hs);

            for(PS::S32 repeat = 0; repeat < 3; repeat++) {
                vndf hi_i  = rcp(hs_i);
                vndf hi3_i = ND::calcVolumeInverse(hi_i);
                vndf hi4_i = hi_i * hi3_i;
                vndf rh_i(0.);
                vndf nj_i(0.);

                for(PS::S32 j = 0; j < njp; j++) {
                    vndf dx_ij = px_i - vndf(epj[j].pos[0]);
                    vndf dy_ij = py_i - vndf(epj[j].pos[1]);
                    vndf dz_ij = pz_i - vndf(epj[j].pos[2]);

                    vndf r2_ij = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;

                    vndf r1_ij = vndf::sqrt(r2_ij);
                    vndf q_i   = r1_ij * hi_i;

                    vndf kw0 = SK::kernel0th(q_i);
                    vndf rhj = vndf(epj[j].mass) * hi3_i * kw0;

                    rh_i += rhj;
#ifdef __AVX512F__
                    nj_i += _mm512_mask_mov_pd(vndf(0.),
                                               _mm512_cmp_pd_mask(q_i, vndf(1.), _CMP_LT_OS),
                                               vndf(1.));
#else
                    nj_i += ((q_i < 1.) & vndf(1.));                    
#endif
                }

                PS::F64 buf_dens[nvector] __attribute__((aligned(vndf::nvector*8)));
                PS::F64 buf_nj[nvector]   __attribute__((aligned(vndf::nvector*8)));
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
            
            vndf hi_i  = rcp(hs_i);
            vndf hi4_i = hi_i * ND::calcVolumeInverse(hi_i);
            vndf grdh_i(0.);
            vndf divv_i(0.);
            vndf rotx_i(0.);
            vndf roty_i(0.);
            vndf rotz_i(0.);
            for(PS::S32 j = 0; j < njp; j++) {
                vndf dpx_ij = px_i - vndf(epj[j].pos[0]);
                vndf dpy_ij = py_i - vndf(epj[j].pos[1]);
                vndf dpz_ij = pz_i - vndf(epj[j].pos[2]);
                vndf dvx_ij = vx_i - vndf(epj[j].vel[0]);
                vndf dvy_ij = vy_i - vndf(epj[j].vel[1]);
                vndf dvz_ij = vz_i - vndf(epj[j].vel[2]);

                vndf r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                vndf ri_ij = rsqrt(r2_ij);
#ifdef __AVX512F__
                ri_ij = _mm512_mask_mov_pd(vndf(0.),
                    _mm512_cmp_pd_mask(id_i, vndf(epj[j].id), _CMP_NEQ_UQ),
                    ri_ij);
#else
                ri_ij = ((id_i != vndf(epj[j].id)) & ri_ij);
#endif
                vndf r1_ij = r2_ij * ri_ij;
                vndf q_i = r1_ij * hi_i;

                vndf kw0 = SK::kernel0th(q_i);
                vndf kw1 = SK::kernel1st(q_i);

                vndf m_j(epj[j].mass);
                vndf ghj = vndf(RP::NumberOfDimension) * kw0;
                ghj += q_i * kw1;
                grdh_i -= ghj * hi4_i * m_j;

                vndf dw_ij  = m_j * hi4_i * kw1 * ri_ij;
                vndf dwx_ij = dw_ij * dpx_ij;
                vndf dwy_ij = dw_ij * dpy_ij;
                vndf dwz_ij = dw_ij * dpz_ij;

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

#ifdef __AVX512F__
            PS::F64 buf_dens[nvector] __attribute__((aligned(vndf::nvector*8)));
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_dens[ii] = density[i+ii].dens;
            }
            vndf dens_i;
            dens_i.load(buf_dens);
#else
            vndf dens_i(density[i].dens, density[i+1].dens, density[i+2].dens, density[i+3].dens);
#endif
            vndf deni_i = rcp(dens_i);
            vndf omgi_i = rcp(vndf(1.) + hs_i * deni_i * grdh_i * rcp(RP::NumberOfDimension));
            vndf rot2_i = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
#ifdef __AVX512F__
            vndf rotv_i = rot2_i * _mm512_mask_mov_pd(vndf(0.),
                _mm512_cmp_pd_mask(rot2_i, vndf(0.), _CMP_NEQ_UQ),
                rsqrt(rot2_i));
#else
            vndf rotv_i = rot2_i * ((rot2_i != 0.) & rsqrt(rot2_i));
#endif
            rotv_i *= deni_i * omgi_i;
            divv_i *= deni_i * omgi_i;

            PS::F64 buf_divv[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_rotv[nvector] __attribute__((aligned(vndf::nvector*8)));
            PS::F64 buf_omgi[nvector] __attribute__((aligned(vndf::nvector*8)));
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

#endif

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
            // tanikawa modified this 171205 FROM
            omgi_i = (hs_i * deni_i * grdh_i / RP::NumberOfDimension != -1.) ?
                omgi_i : 1.;
            // tanikawa modified this 171205 TO
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
