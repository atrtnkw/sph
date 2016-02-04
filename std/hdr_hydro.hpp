#pragma once

class Hydro {
public:
    PS::F64vec acch;
    PS::F64vec accg;
    PS::F64    udot;
    PS::F64    vsmx;
    PS::F64    diffu;
    void clear() {
        this->acch  = 0.;
        this->accg  = 0.;
        this->udot  = 0.;
        this->vsmx  = 0.;
        this->diffu = 0.;
    }
};

void SPH::copyFromForce(const Hydro & hydro) {
    using namespace CodeUnit;
    this->acch  = hydro.acch;
    this->accg1 = hydro.accg * GravityConstantInThisUnit;
    this->udot  = hydro.udot;
    this->vsmx  = hydro.vsmx;
    this->diffu = hydro.diffu;
}

class HydroEPI {
public:
    PS::S64    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    ksr;
    PS::F64    hinv;
    PS::F64    pres;
    PS::F64    presc;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;
    PS::F64    alphu;
    PS::F64    thrm;
    PS::F64    eta;

    void copyFromFP(const SPH & sph){ 
        this->id    = sph.id;
        this->pos   = sph.pos;
        this->vel   = sph.vel;
        this->ksr   = sph.ksr;
        this->hinv  = 1. / sph.ksr;
        this->pres  = sph.pres;
        this->presc = sph.pres * sph.grdh / (sph.dens * sph.dens);
        this->bswt  = 0.5 * sph.bswt;
        this->dens  = sph.dens;
        this->vsnd  = sph.vsnd;
        this->alph  = sph.alph;
        this->alphu = sph.alphu;        
        //PS::F64 ene = sph.uene - CalcEquationOfState::getEnergyMin(sph.dens, sph.abar, sph.zbar);
        PS::F64 ene = sph.uene - sph.umin;
        this->thrm  = std::max(ene, 0.);
        this->eta   = sph.eta;
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
    PS::F64    pres;
    PS::F64    presc;
    PS::F64    bswt;
    PS::F64    dens;
    PS::F64    vsnd;
    PS::F64    alph;
    PS::F64    alphu;
    PS::F64    thrm;
    PS::F64    eta;

    void copyFromFP(const SPH & sph){ 
        this->id    = sph.id;
        this->mass  = sph.mass;
        this->pos   = sph.pos;
        this->vel   = sph.vel;
        this->ksr   = sph.ksr;
        this->hinv  = 1. / sph.ksr;
        this->pres  = sph.pres;
        this->presc = sph.pres * sph.grdh / (sph.dens * sph.dens);
        this->bswt  = 0.5 * sph.bswt;
        this->dens  = sph.dens;
        this->vsnd  = sph.vsnd;
        this->alph  = sph.alph;
        this->alphu = sph.alphu;
        //PS::F64 ene = sph.uene - CalcEquationOfState::getEnergyMin(sph.dens, sph.abar, sph.zbar);
        PS::F64 ene = sph.uene - sph.umin;
        this->thrm  = std::max(ene, 0.);
        this->eta   = sph.eta;
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
            PS::F64 buf_pp[nvector];
            PS::F64 buf_pc[nvector];
            PS::F64 buf_bs[nvector];
            PS::F64 buf_dn[nvector];
            PS::F64 buf_cs[nvector];
            PS::F64 buf_al[nvector];
            PS::F64 buf_au[nvector];
            PS::F64 buf_eg[nvector];
            PS::F64 buf_et[nvector];
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
                buf_pp[ii] = epi[i+ii].pres;
                buf_pc[ii] = epi[i+ii].presc;
                buf_bs[ii] = epi[i+ii].bswt;
                buf_dn[ii] = epi[i+ii].dens;
                buf_cs[ii] = epi[i+ii].vsnd;
                buf_al[ii] = epi[i+ii].alph;
                buf_au[ii] = epi[i+ii].alphu;
                buf_eg[ii] = epi[i+ii].thrm;
                buf_et[ii] = epi[i+ii].eta;
            }
            v4df id_i;
            v4df px_i, py_i, pz_i;
            v4df vx_i, vy_i, vz_i;
            v4df hi_i, pp_i, pc_i;
            v4df bs_i, dn_i, cs_i, al_i, au_i, eg_i;
            v4df et_i;
            id_i.load(buf_id);
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            vx_i.load(buf_vx);
            vy_i.load(buf_vy);
            vz_i.load(buf_vz);
            hi_i.load(buf_hi);
            pp_i.load(buf_pp);
            pc_i.load(buf_pc);
            bs_i.load(buf_bs);
            dn_i.load(buf_dn);
            cs_i.load(buf_cs);
            al_i.load(buf_al);
            au_i.load(buf_au);
            eg_i.load(buf_eg);
            et_i.load(buf_et);
            v4df hi4_i = hi_i * ND::calcVolumeInverse(hi_i);
            v4df achx_i(0.);
            v4df achy_i(0.);
            v4df achz_i(0.);
            v4df udot_i(0.);
            v4df vsmx_i(0.);
            v4df difu_i(0.);
            v4df g1x_i(0.);
            v4df g1y_i(0.);
            v4df g1z_i(0.);
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
                v4df dw_i  = hi4_i * SK::kernel1st(q_i);
                v4df dw_j  = hi4_j * SK::kernel1st(q_j);
                v4df ka_ij = (dw_i + dw_j) * v4df(epj[j].mass);

                v4df rv_ij  = dpx_ij * dvx_ij + dpy_ij * dvy_ij + dpz_ij * dvz_ij;
                v4df w_ij   = rv_ij * ri_ij;
                v4df w0_ij  = ((w_ij < v4df(0.)) & w_ij);
                v4df vs_ij  = cs_i + v4df(epj[j].vsnd) - v4df(3.) * w0_ij;
                v4df rhi_ij = rcp(dn_i + v4df(epj[j].dens));
                v4df av0_ij = (bs_i + v4df(epj[j].bswt)) * (al_i + v4df(epj[j].alph))
                    * vs_ij * w0_ij;
                                
                vsmx_i  = v4df::max(vsmx_i, vs_ij);

                v4df ta_ij = (pc_i + v4df(epj[j].presc) - v4df(0.5) * av0_ij * rhi_ij)
                    * ka_ij * ri_ij;
                achx_i -= ta_ij * dpx_ij;
                achy_i -= ta_ij * dpy_ij;
                achz_i -= ta_ij * dpz_ij;

                v4df vsu2_ij  = v4df::fabs(pp_i - v4df(epj[j].pres)) * rhi_ij * v4df(2.);
                v4df vsui_ij  = ((vsu2_ij != 0.) & rsqrt(vsu2_ij));
                v4df vsu_ij   = vsu2_ij * vsui_ij;
                v4df du_ij    = eg_i - v4df(epj[j].thrm);

                udot_i += ka_ij * (pc_i * w_ij
                                   - rhi_ij * (v4df(0.25) * av0_ij * w0_ij
                                               - (au_i + v4df(epj[j].alphu)) * vsu_ij * du_ij));

                difu_i += rhi_ij * du_ij * ka_ij * ri_ij;

                v4df dg_ij = v4df(epj[j].mass) * ri_ij * (et_i * dw_i + v4df(epj[j].eta) * dw_j);
                g1x_i += dg_ij * dpx_ij;
                g1y_i += dg_ij * dpy_ij;
                g1z_i += dg_ij * dpz_ij;
                    
            }

            achx_i *= v4df(0.5);
            achy_i *= v4df(0.5);
            achz_i *= v4df(0.5);
            udot_i *= v4df(0.5);
            difu_i *= v4df(2.);
            g1x_i  *= v4df(0.5);
            g1y_i  *= v4df(0.5);
            g1z_i  *= v4df(0.5);

            PS::F64 buf_ax[nvector];
            PS::F64 buf_ay[nvector];
            PS::F64 buf_az[nvector];
            PS::F64 buf_du[nvector];
            PS::F64 buf_vs[nvector];
            PS::F64 buf_df[nvector];
            PS::F64 buf_gx[nvector];
            PS::F64 buf_gy[nvector];
            PS::F64 buf_gz[nvector];
            achx_i.store(buf_ax);
            achy_i.store(buf_ay);
            achz_i.store(buf_az);
            udot_i.store(buf_du);
            vsmx_i.store(buf_vs);
            difu_i.store(buf_df);
            g1x_i.store(buf_gx);
            g1y_i.store(buf_gy);
            g1z_i.store(buf_gz);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                hydro[i+ii].acch[0] = buf_ax[ii];
                hydro[i+ii].acch[1] = buf_ay[ii];
                hydro[i+ii].acch[2] = buf_az[ii];
                hydro[i+ii].udot    = buf_du[ii];
                hydro[i+ii].vsmx    = buf_vs[ii];
                hydro[i+ii].diffu   = buf_df[ii];
                hydro[i+ii].accg[0] = buf_gx[ii];
                hydro[i+ii].accg[1] = buf_gy[ii];
                hydro[i+ii].accg[2] = buf_gz[ii];
            }
            
        }

    }
};
