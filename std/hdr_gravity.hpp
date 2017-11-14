namespace ParticleSimulator {

    class GravityMonopole{
    public:
        F32    mass;
        F32vec pos;
        F32    eps2;
        GravityMonopole(){
            mass = 0.0;
            pos  = 0.0;
            eps2 = 0.0;
        }
        GravityMonopole(const F32 m,
                        const F32vec & p,
                        const F32 e2){
            mass = m;
            pos  = p;
            eps2 = e2;
        }
        void init(){
            mass = 0.0;
            pos  = 0.0;
            eps2 = 0.0;
        }
        F32vec getPos() const {
            return pos;
        }
        F32 getCharge() const {
            return mass;
        }
        F32 getSofteningSquare() const {
            return eps2;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            mass += epj.getCharge();
            pos  += epj.getCharge() * epj.getPos();
            eps2 += epj.getCharge() * epj.getSofteningSquare();
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            F32 minv = 1. / mass;
            pos  = pos  * minv;
            eps2 = eps2 * minv;
        }
        void accumulate(const GravityMonopole & mom){
            mass += mom.getCharge();
            pos  += mom.getCharge() * mom.getPos();
            eps2 += mom.getCharge() * mom.getSofteningSquare();
        }
        void accumulate2(const GravityMonopole & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class GravitySPJ {
    public:
        F64    mass;
        F64vec pos;
        F64    eps2;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            mass = mom.getCharge();
            pos  = mom.getPos();
            eps2 = mom.getSofteningSquare();
        }
        void clear(){
            mass = 0.0;
            pos  = 0.0;
            eps2 = 0.0;
        }
        F64 getCharge() const {
            return mass;
        }
        F64vec getPos() const {
            return pos;
        }
        void setPos(const F64vec & pos_new) {
            pos = pos_new;
        }
        GravityMonopole convertToMoment() const {
            return GravityMonopole(mass, pos, eps2);
        }
    };

}

class Gravity{
public:
    PS::F64vec acc;
    PS::F64    pot;
    PS::F64    eta;
    void clear(){
        acc = 0.0;
        pot = 0.0;
        eta = 0.0;
    }
};

void SPH::copyFromForce(const Gravity & gravity) {
    using namespace CodeUnit;
    this->accg2 = GravityConstantInThisUnit * gravity.acc;
    this->eta   = this->ksr * this->ksr * SK::ksrhinv * SK::ksrhinv
        * this->grdh / (SK::dim * this->dens) * gravity.eta;
    this->pot   = GravityConstantInThisUnit * gravity.pot;
}

void MassLess::copyFromForce(const Gravity & gravity) {
    using namespace CodeUnit;
    this->accg2 = GravityConstantInThisUnit * gravity.acc;
    this->eta   = this->ksr * this->ksr * SK::ksrhinv * SK::ksrhinv
        * this->grdh / (SK::dim * this->dens) * gravity.eta;
    this->pot   = GravityConstantInThisUnit * gravity.pot;
}

void BlackHoleNeutronStar::copyFromForce(const Gravity & gravity) {
    using namespace CodeUnit;
    this->acc   = GravityConstantInThisUnit * gravity.acc;
    this->eta   = 0.;
    this->pot   = GravityConstantInThisUnit * gravity.pot;
}

class GravityEPI {
public:
    PS::F64vec pos;
    PS::F64    eps2;
    void copyFromFP(const SPH & sph){ 
        pos  = sph.pos;
        eps2 = sph.ksr * sph.ksr * SK::ksrhinv * SK::ksrhinv;
    }

    void copyFromFP(const MassLess & msls){         
        pos  = msls.pos;
        eps2 = RP::GravitationalSoftening * RP::GravitationalSoftening;
    }

    void copyFromFP(const BlackHoleNeutronStar & bhns) {
        pos  = bhns.pos;
        eps2 = bhns.eps * bhns.eps;
    }

    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
};

class GravityEPJ {
public:
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    eps2;
    void copyFromFP(const SPH & sph){ 
        mass = sph.mass;
        pos  = sph.pos;
        eps2 = sph.ksr * sph.ksr * SK::ksrhinv * SK::ksrhinv;
    }

    void copyFromFP(const MassLess & msls){ 
        mass = msls.mass;
        pos  = msls.pos;
        eps2 = RP::GravitationalSoftening * RP::GravitationalSoftening;
    }

    void copyFromFP(const BlackHoleNeutronStar & bhns) {
        mass = bhns.mass;
        pos  = bhns.pos;
        eps2 = bhns.eps * bhns.eps;
    }

    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64 getSofteningSquare() const {
        return this->eps2;
    }
};

#ifdef USE_INTRINSICS

template <class TGravityJ>
struct calcGravity {    

    void operator () (const GravityEPI *epi,
                      const PS::S32 nip,
                      const TGravityJ *epj,
                      const PS::S32 njp,
                      Gravity *gravity) {

        v4df (*rsqrt)(v4df) = v4df::rsqrt_4th;
        PS::S32 nvector = v4df::getVectorLength();

        for(PS::S32 i = 0; i < nip; i += nvector) {
            PS::F64 buf_px[nvector] __attribute__((aligned(32)));
            PS::F64 buf_py[nvector] __attribute__((aligned(32)));
            PS::F64 buf_pz[nvector] __attribute__((aligned(32)));
            PS::F64 buf_e2[nvector] __attribute__((aligned(32)));
            const PS::S32 nii = std::min(nip - i, nvector);
            for(PS::S32 ii = 0; ii < nii; ii++) {
                buf_px[ii] = epi[i+ii].pos[0];
                buf_py[ii] = epi[i+ii].pos[1];
                buf_pz[ii] = epi[i+ii].pos[2];
                buf_e2[ii] = epi[i+ii].eps2;
            }
            v4df px_i;
            v4df py_i;
            v4df pz_i;
            v4df e2_i;
            px_i.load(buf_px);
            py_i.load(buf_py);
            pz_i.load(buf_pz);
            e2_i.load(buf_e2);
            v4df ax_i(0.0);
            v4df ay_i(0.0);
            v4df az_i(0.0);
            v4df pt_i(0.0);
            v4df et_i(0.0);
            for(PS::S32 j = 0; j < njp; j++) {
                v4df dpx_ij = px_i - v4df(epj[j].pos[0]);
                v4df dpy_ij = py_i - v4df(epj[j].pos[1]);
                v4df dpz_ij = pz_i - v4df(epj[j].pos[2]);

                v4df r2_ij  = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                v4df re2_i  = r2_ij + e2_i;
                v4df rei_i  = rsqrt(re2_i);
                v4df rei3_i = rei_i * rei_i * rei_i;
                v4df re2_j  = r2_ij + v4df(epj[j].eps2);
                v4df rei_j  = rsqrt(re2_j);
                v4df rei3_j = rei_j * rei_j * rei_j;

                v4df m_j    = v4df(epj[j].mass);
                v4df dg2_ij = m_j * (rei3_i + rei3_j);

                pt_i -= m_j    * (rei_i  + rei_j);
                ax_i -= dpx_ij * dg2_ij;
                ay_i -= dpy_ij * dg2_ij;
                az_i -= dpz_ij * dg2_ij;
                et_i += m_j    * rei3_i;

            }
            ax_i *= v4df(0.5);
            ay_i *= v4df(0.5);
            az_i *= v4df(0.5);
            pt_i *= v4df(0.5);

            PS::F64 buf_ax[nvector] __attribute__((aligned(32)));
            PS::F64 buf_ay[nvector] __attribute__((aligned(32)));
            PS::F64 buf_az[nvector] __attribute__((aligned(32)));
            PS::F64 buf_pt[nvector] __attribute__((aligned(32)));
            PS::F64 buf_et[nvector] __attribute__((aligned(32)));

            ax_i.store(buf_ax);
            ay_i.store(buf_ay);
            az_i.store(buf_az);
            pt_i.store(buf_pt);
            et_i.store(buf_et);

            for(PS::S32 ii = 0; ii < nii; ii++) {
                gravity[i+ii].acc[0] += buf_ax[ii];
                gravity[i+ii].acc[1] += buf_ay[ii];
                gravity[i+ii].acc[2] += buf_az[ii];
                gravity[i+ii].pot    += buf_pt[ii];
                gravity[i+ii].eta    += buf_et[ii];
            }
        }
    }
};

#else

template <class TGravityJ>
struct calcGravity {    

    void operator () (const GravityEPI *epi,
                      const PS::S32 nip,
                      const TGravityJ *epj,
                      const PS::S32 njp,
                      Gravity *gravity) {

        for(PS::S32 i = 0; i < nip; i++) {
            PS::F64 px_i = epi[i].pos[0];
            PS::F64 py_i = epi[i].pos[1];
            PS::F64 pz_i = epi[i].pos[2];
            PS::F64 e2_i = epi[i].eps2;
            PS::F64 ax_i = 0.;
            PS::F64 ay_i = 0.;
            PS::F64 az_i = 0.;
            PS::F64 pt_i = 0.;
            PS::F64 et_i = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::F64 dpx_ij = px_i - epj[j].pos[0];
                PS::F64 dpy_ij = py_i - epj[j].pos[1];
                PS::F64 dpz_ij = pz_i - epj[j].pos[2];

                PS::F64 r2_ij  = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
                PS::F64 re2_i  = r2_ij + e2_i;
                PS::F64 rei_i  = 1. / sqrt(re2_i);
                PS::F64 rei3_i = rei_i * rei_i * rei_i;
                PS::F64 re2_j  = r2_ij + epj[j].eps2;
                PS::F64 rei_j  = 1. / sqrt(re2_j);
                PS::F64 rei3_j = rei_j * rei_j * rei_j;

                PS::F64 m_j    = epj[j].mass;
                PS::F64 dg2_ij = m_j * (rei3_i + rei3_j);
                
                pt_i -= m_j    * (rei_i  + rei_j);
                ax_i -= dpx_ij * dg2_ij;
                ay_i -= dpy_ij * dg2_ij;
                az_i -= dpz_ij * dg2_ij;
                et_i += m_j    * rei3_i;
            }

            ax_i *= 0.5;
            ay_i *= 0.5;
            az_i *= 0.5;
            pt_i *= 0.5;

            gravity[i].acc[0] += ax_i;
            gravity[i].acc[1] += ay_i;
            gravity[i].acc[2] += az_i;
            gravity[i].pot    += pt_i;
            gravity[i].eta    += et_i;

        }

    }
};

#endif

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls,
          class Tgravity>
void calcGravityKernel(Tdinfo & dinfo,
                       Tsph & sph,
                       Tbhns & bhns,
                       Tmsls & msls,
                       Tgravity & gravity) {
    if(RP::FlagGravity == 0) {
        return;
    }

    gravity.setParticleLocalTree(sph);
    gravity.setParticleLocalTree(msls, false);
    gravity.setParticleLocalTree(bhns, false);
    gravity.calcForceMakingTree(calcGravity<GravityEPJ>(),
                                calcGravity<PS::GravitySPJ>(),
                                dinfo);
    PS::S32 nsph  = sph.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nsph; i++) {
        sph[i].copyFromForce(gravity.getForce(i));
    }
    PS::S32 nmsls = msls.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nmsls; i++) {
        msls[i].copyFromForce(gravity.getForce(i+nsph));
    }
    PS::S32 nbhns = bhns.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nbhns; i++) {
        bhns[i].copyFromForce(gravity.getForce(i+nsph+nmsls));
    }
}
