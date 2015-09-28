#pragma once

#include <limits>

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

struct calcDensityBasic {

    void operator () (const DensityEPI *epi,
                      const PS::S32 nip,
                      const DensityEPJ *epj,
                      const PS::S32 njp,
                      Density *density) {

        for(PS::S32 i = 0; i < nip; i++) {
            PS::S32    id_i = epi[i].id;
            PS::F64vec x_i  = epi[i].pos;
            PS::F64vec v_i  = epi[i].vel;
            PS::F64    h_i  = epi[i].ksr;

            for(PS::S32 repeat = 0; repeat < 3; repeat++) {
                PS::F64 hi_i  = 1.d / h_i;
                PS::F64 hi3_i = SPH::calcVolumeInverse(hi_i);
                PS::F64 hi4_i = hi_i * hi3_i;
                PS::F64 rh_i  = 0.;
                PS::S32 nj_i  = 0;

                for(PS::S32 j = 0; j < njp; j++) {
                    PS::F64    m_j = epj[j].mass;
                    PS::F64vec x_j = epj[j].pos;
                    
                    PS::F64vec dx_ij = x_i - x_j;
                    PS::F64    r2_ij = dx_ij * dx_ij;

                    PS::F64    r1_ij = sqrt(r2_ij);
                    PS::F64    q_i   = r1_ij * hi_i;

                    PS::F64 kw0 = KernelSph::kernel0th(q_i);
                    
                    PS::F64 rhj =   m_j * hi3_i * kw0;
                    
                    rh_i   += rhj;
                    nj_i   += (q_i < 1.d) ? 1 : 0;

                }
                density[i].dens = rh_i;
                density[i].np   = nj_i;

                h_i = KernelSph::eta * KernelSph::ksrh
                    * SPH::calcPowerOfDimInverse(epi[i].mass, rh_i);
                density[i].ksr = h_i;
                density[i].itr = (h_i > epi[i].rs) ? true : false;
            }

            PS::F64    hi_i   = 1.d / h_i;
            PS::F64    hi4_i  = hi_i * SPH::calcVolumeInverse(hi_i);
            PS::F64    gh_i  = 0.;
            PS::F64    divv_i = 0.;
            PS::F64vec rotv_i = 0.;            

            PS::F64    h2c2_i = h_i * h_i * KernelSph::ksrhinv * KernelSph::ksrhinv;

            for(PS::S32 j = 0; j < njp; j++) {
                PS::S32    id_j = epj[j].id;
                PS::F64    m_j  = epj[j].mass;
                PS::F64vec x_j  = epj[j].pos;
                PS::F64vec v_j  = epj[j].vel;

                PS::F64vec dx_ij = x_i - x_j;
                PS::F64vec dv_ij = v_i - v_j;
                PS::F64    r2_ij = dx_ij * dx_ij;
                PS::F64    r1_ij = sqrt(r2_ij);
                PS::F64    ri_ij = (id_i != id_j) ? 1. / r1_ij : 0.;
                PS::F64    q_i   = r1_ij * hi_i;

                PS::F64 kw0 = KernelSph::kernel0th(q_i);
                PS::F64 kw1 = KernelSph::kernel1st(q_i);
                
                PS::F64    ghj   = - m_j * hi4_i * (KernelSph::dim * kw0 + q_i * kw1);
                PS::F64vec dw_ij = (m_j * hi4_i * kw1 * ri_ij) * dx_ij;

                gh_i   += ghj;
                divv_i -= dv_ij * dw_ij;
                rotv_i += dv_ij ^ dw_ij;

            }

            PS::F64 rhi_i = 1. / density[i].dens;
            density[i].grdh = 1.d / (1.d + h_i * rhi_i * gh_i / KernelSph::dim);
            PS::F64 grd_i = density[i].grdh;
            density[i].rotv = sqrt(rotv_i * rotv_i) * rhi_i * grd_i;
            density[i].divv = divv_i * rhi_i * grd_i;
        }        
    }
};

struct calcDensityX86 {

    void operator () (const DensityEPI *epi,
                      const PS::S32 nip,
                      const DensityEPJ *epj,
                      const PS::S32 njp,
                      Density *density) {

#ifdef TIMETEST
        ninteract += nip * njp;
#endif

        v4df (*rcp)(v4df)   = v4df::rcp_4th;
        v4df (*rsqrt)(v4df) = v4df::rsqrt_4th;

        const PS::S32 nvector = v4df::getVectorLength();
        
        for(PS::S32 i = 0; i < nip; i += nvector) {

            const PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;

            PS::F64 buf_id[nvector];
            PS::F64 buf_px[nvector];
            PS::F64 buf_py[nvector];
            PS::F64 buf_pz[nvector];
            PS::F64 buf_vx[nvector];
            PS::F64 buf_vy[nvector];
            PS::F64 buf_vz[nvector];
            PS::F64 buf_h[nvector];
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_id[ii] = epi[i+ii].id;
                buf_px[ii] = epi[i+ii].pos[0];
                buf_py[ii] = epi[i+ii].pos[1];
                buf_pz[ii] = epi[i+ii].pos[2];
                buf_vx[ii] = epi[i+ii].vel[0];
                buf_vy[ii] = epi[i+ii].vel[1];
                buf_vz[ii] = epi[i+ii].vel[2];
                buf_h[ii]  = epi[i+ii].ksr;
            }
            v4df id_i;
            v4df px_i;
            v4df py_i;
            v4df pz_i;
            v4df vx_i;
            v4df vy_i;
            v4df vz_i;
            v4df h_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            h_i.load(buf_h);

            for(PS::S32 repeat = 0; repeat < 3; repeat++) {
                v4df hi_i  = rcp(h_i);
                v4df hi3_i = SPH::calcVolumeInverse(hi_i);
                v4df hi4_i = hi_i * hi3_i;
                v4df rh_i(0.d);
                v4df nj_i(0.d);

#ifdef TIMETEST
                PS::F64 t1 = getWallclockTime();;
#endif
                for(PS::S32 j = 0; j < njp; j++) {
                    v4df dx_ij = px_i - v4df(epj[j].pos[0]);
                    v4df dy_ij = py_i - v4df(epj[j].pos[1]);
                    v4df dz_ij = pz_i - v4df(epj[j].pos[2]);

                    v4df r2_ij = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;

#if 0
                    v4df r1_ij = r2_ij * rsqrt(r2_ij);
                    //r1_ij = ((id_i != v4df(epj[j].id)) & r1_ij);
                    r1_ij = ((r2_ij != v4df(0.d)) & r1_ij);
#else
                    v4df r1_ij = v4df::sqrt(r2_ij);
#endif
                    v4df q_i   = r1_ij * hi_i;

                    v4df kw0 = KernelSph::kernel0th(q_i);
                    v4df rhj = v4df(epj[j].mass) * hi3_i * kw0;

                    rh_i += rhj;
                    nj_i += ((q_i < 1.d) & v4df(1.d));                    
                }
#ifdef TIMETEST
                tcalcdens += getWallclockTime() - t1;
                ncalcdens++;
#endif

                PS::F64 buf0[nvector], buf1[nvector], hs[nvector];
                rh_i.store(buf0);
                nj_i.store(buf1);
                for(PS::S32 ii = 0; ii < nii; ii++) {
                    hs[ii] = KernelSph::eta * KernelSph::ksrh
                        * SPH::calcPowerOfDimInverse(epi[i+ii].mass, buf0[ii]);
                    hs[ii] = std::min(hs[ii], SPH::ksrmax);
                    density[i+ii].dens = buf0[ii];
                    density[i+ii].np   = (PS::S32)buf1[ii];
                    density[i+ii].ksr  = hs[ii];
                    density[i+ii].itr  = (hs[ii] > epi[i+ii].rs) ? true : false;
                }
                h_i.load(hs);
            }

            h_i  = v4df(density[i].ksr, density[i+1].ksr, density[i+2].ksr, density[i+3].ksr);
            v4df hi_i  = rcp(h_i);
            v4df hi4_i = hi_i * SPH::calcVolumeInverse(hi_i);
            v4df gh_i(0.d);
            v4df divv_i(0.d);
            v4df rotx_i(0.d);            
            v4df roty_i(0.d);
            v4df rotz_i(0.d);

#ifdef TIMETEST
            PS::F64 t1 = getWallclockTime();;
#endif
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

                v4df kw0 = KernelSph::kernel0th(q_i);
                v4df kw1 = KernelSph::kernel1st(q_i);

                v4df m_j(epj[j].mass);
                v4df ghj = v4df(KernelSph::dim) * kw0;
                ghj += q_i * kw1;
                gh_i -= ghj * hi4_i * m_j;

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
#ifdef TIMETEST
            tcalcgrdh += getWallclockTime() - t1;
            ncalcgrdh++;
#endif

            v4df rh_i(density[i].dens, density[i+1].dens, density[i+2].dens, density[i+3].dens);
            v4df rhi_i = rcp(rh_i);
            v4df grd_i = rcp(v4df(1.d) + h_i * rhi_i * gh_i * rcp(KernelSph::dim));
            v4df rot2  = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
            v4df roti  = rsqrt(rot2);
            v4df rot1  = rot2 * ((rot2 != 0.d) & roti);
            v4df rotv  = rot1 * rhi_i * grd_i;
            v4df divv  = divv_i * rhi_i * grd_i;

            PS::F64 buf0[nvector], buf1[nvector], buf2[nvector];
            grd_i.store(buf0);
            rotv.store(buf1);
            divv.store(buf2);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                density[i+ii].grdh = buf0[ii];
                density[i+ii].rotv = buf1[ii];
                density[i+ii].divv = buf2[ii];
            }

        }        
    }
};

#if 0
struct calcDensitySingleX86 {

    void operator () (const DensityEPI *epi,
                      const PS::S32 nip,
                      const DensityEPJ *epj,
                      const PS::S32 njp,
                      Density *density) {

#ifdef TIMETEST
        ninteract += nip * njp;
#endif

        v8sf (*rcp)(v8sf)   = v8sf::rcp_1st;
        v8sf (*rsqrt)(v8sf) = v8sf::rsqrt_1st;

        const PS::S32 nvector = v8sf::getVectorLength();
        
        for(PS::S32 i = 0; i < nip; i += nvector) {

            const PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;

            v8sf id_i(epi[i].id,     epi[i+1].id,     epi[i+2].id,     epi[i+3].id,
                      epi[i+4].id,   epi[i+5].id,     epi[i+6].id,     epi[i+7].id);
            v8sf px_i(epi[i].pos[0],   epi[i+1].pos[0], epi[i+2].pos[0], epi[i+3].pos[0],
                      epi[i+4].pos[0], epi[i+5].pos[0], epi[i+6].pos[0], epi[i+7].pos[0]);
            v8sf py_i(epi[i].pos[1],   epi[i+1].pos[1], epi[i+2].pos[1], epi[i+3].pos[1],
                      epi[i+4].pos[1], epi[i+5].pos[1], epi[i+6].pos[1], epi[i+7].pos[1]);
            v8sf pz_i(epi[i].pos[2],   epi[i+1].pos[2], epi[i+2].pos[2], epi[i+3].pos[2],
                      epi[i+4].pos[2], epi[i+5].pos[2], epi[i+6].pos[2], epi[i+7].pos[2]);
            v8sf vx_i(epi[i].vel[0],   epi[i+1].vel[0], epi[i+2].vel[0], epi[i+3].vel[0],
                      epi[i+4].vel[0], epi[i+5].vel[0], epi[i+6].vel[0], epi[i+7].vel[0]);
            v8sf vy_i(epi[i].vel[1],   epi[i+1].vel[1], epi[i+2].vel[1], epi[i+3].vel[1],
                      epi[i+4].vel[1], epi[i+5].vel[1], epi[i+6].vel[1], epi[i+7].vel[1]);
            v8sf vz_i(epi[i].vel[2],   epi[i+1].vel[2], epi[i+2].vel[2], epi[i+3].vel[2],
                      epi[i+4].vel[2], epi[i+5].vel[2], epi[i+6].vel[2], epi[i+7].vel[2]);
            v8sf h_i(epi[i].ksr,   epi[i+1].ksr, epi[i+2].ksr, epi[i+3].ksr,
                     epi[i+4].ksr, epi[i+5].ksr, epi[i+6].ksr, epi[i+7].ksr);

            for(PS::S32 repeat = 0; repeat < 3; repeat++) {
                v8sf hi_i  = rcp(h_i);
                v8sf hi3_i = SPH::calcVolumeInverse(hi_i);
                v8sf hi4_i = hi_i * hi3_i;
                v8sf rh_i(0.d);
                v8sf nj_i(0.d);

#ifdef TIMETEST
                PS::F64 t1 = getWallclockTime();;
#endif
                for(PS::S32 j = 0; j < njp; j++) {
                    v8sf dx_ij = px_i - v8sf(epj[j].pos[0]);
                    v8sf dy_ij = py_i - v8sf(epj[j].pos[1]);
                    v8sf dz_ij = pz_i - v8sf(epj[j].pos[2]);

                    v8sf r2_ij = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;

                    v8sf r1_ij = r2_ij * rsqrt(r2_ij);
                    r1_ij = ((id_i != v8sf(epj[j].id)) & r1_ij);
                    v8sf q_i   = r1_ij * hi_i;

                    v8sf kw0 = KernelSph::kernel0th(q_i);
                    v8sf rhj = v8sf(epj[j].mass) * hi3_i * kw0;

                    rh_i += rhj;
                    nj_i += ((q_i < 1.d) & v8sf(1.d));
                    
                    // vcvt       x 3
                    // vbroadcast x 5
                    // vmova      x 1
                    // vcmp       x 2
                    // vand       x 2
                    // vmax       x 2
                    // vrsqrt     x 1
                    // vadd       x 6
                    // vmul       x 11
                    // vfmadd     x 9 
                }
#ifdef TIMETEST
                tcalcdens += getWallclockTime() - t1;
                ncalcdens++;
#endif

                PS::F32 buf0[nvector], buf1[nvector], hs[nvector];
                rh_i.store(buf0);
                nj_i.store(buf1);
                for(PS::S32 ii = 0; ii < nii; ii++) {
                    hs[ii] = KernelSph::eta * KernelSph::ksrh
                        * SPH::calcPowerOfDimInverse(epi[i+ii].mass, buf0[ii]);
                    density[i+ii].dens = buf0[ii];
                    density[i+ii].np   = (PS::S32)buf1[ii];
                    density[i+ii].ksr  = hs[ii];
                    density[i+ii].itr  = (hs[ii] > epi[i+ii].rs) ? true : false;
                }
                h_i.load(hs);
            }

            h_i  = v8sf(density[i].ksr,   density[i+1].ksr, density[i+2].ksr, density[i+3].ksr,
                        density[i+4].ksr, density[i+5].ksr, density[i+6].ksr, density[i+7].ksr);
            v8sf hi_i  = rcp(h_i);
            v8sf hi4_i = hi_i * SPH::calcVolumeInverse(hi_i);
            v8sf gh_i(0.d);
            v8sf divv_i(0.d);
            v8sf rotx_i(0.d);
            v8sf roty_i(0.d);
            v8sf rotz_i(0.d);

#ifdef TIMETEST
            PS::F64 t1 = getWallclockTime();;
#endif
            for(PS::S32 j = 0; j < njp; j++) {
                v8sf dpx_ij = px_i - v8sf(epj[j].pos[0]);
                v8sf dpy_ij = py_i - v8sf(epj[j].pos[1]);
                v8sf dpz_ij = pz_i - v8sf(epj[j].pos[2]);
                v8sf dvx_ij = vx_i - v8sf(epj[j].vel[0]);
                v8sf dvy_ij = vy_i - v8sf(epj[j].vel[1]);
                v8sf dvz_ij = vz_i - v8sf(epj[j].vel[2]);

                v8sf r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v8sf ri_ij = rsqrt(r2_ij);
                ri_ij = ((id_i != v8sf(epj[j].id)) & ri_ij);
                v8sf r1_ij = r2_ij * ri_ij;
                v8sf q_i = r1_ij * hi_i;

                v8sf kw0 = KernelSph::kernel0th(q_i);
                v8sf kw1 = KernelSph::kernel1st(q_i);

                v8sf m_j(epj[j].mass);
                v8sf ghj = v8sf(KernelSph::dim) * kw0;
                ghj += q_i * kw1;
                gh_i -= ghj * hi4_i * m_j;

                v8sf dw_ij  = m_j * hi4_i * kw1 * ri_ij;
                v8sf dwx_ij = dw_ij * dpx_ij;
                v8sf dwy_ij = dw_ij * dpy_ij;
                v8sf dwz_ij = dw_ij * dpz_ij;

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
#ifdef TIMETEST
            tcalcgrdh += getWallclockTime() - t1;
            ncalcgrdh++;
#endif

            v8sf rh_i(density[i].dens,   density[i+1].dens, density[i+2].dens, density[i+3].dens,
                      density[i+4].dens, density[i+5].dens, density[i+6].dens, density[i+7].dens);
            v8sf rhi_i = rcp(rh_i);
            v8sf grd_i = rcp(v8sf(1.d) + h_i * rhi_i * gh_i * rcp(KernelSph::dim));
            v8sf rot2  = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
            v8sf roti  = rsqrt(rot2);
            v8sf rot1  = rot2 * ((rot2 != 0.d) & roti);
            v8sf rotv  = rot1 * rhi_i * grd_i;
            v8sf divv  = divv_i * rhi_i * grd_i;

            PS::F32 buf0[nvector], buf1[nvector], buf2[nvector];
            grd_i.store(buf0);
            rotv.store(buf1);
            divv.store(buf2);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                density[i+ii].grdh = buf0[ii];
                density[i+ii].rotv = buf1[ii];
                density[i+ii].divv = buf2[ii];
            }

        }        
    }
};
#endif

#ifdef ENABLE_SIMDX86
typedef calcDensityX86 calcDensity;
#elif defined ENABLE_SIMDX86_SINGLE
typedef calcDensitySingleX86 calcDensity;
#else
typedef calcDensityBasic calcDensity;
#endif
