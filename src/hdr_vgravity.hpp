namespace ParticleSimulator {

    class MomentMonopoleSymmetrized{
    public:
        F32    mass;
        F32vec pos;
        F32    eps2;
        MomentMonopoleSymmetrized(){
            mass = 0.0;
            pos  = 0.0;
            eps2 = 0.0;
        }
        MomentMonopoleSymmetrized(const F32 m,
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
        void accumulate(const MomentMonopoleSymmetrized & mom){
            mass += mom.getCharge();
            pos  += mom.getCharge() * mom.getPos();
            eps2 += mom.getCharge() * mom.getSofteningSquare();
        }
        void accumulate2(const MomentMonopoleSymmetrized & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class SPJMonopoleSymmetrized {
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
        MomentMonopoleSymmetrized convertToMoment() const {
            return MomentMonopoleSymmetrized(mass, pos, eps2);
        }
    };

}

class SymmetrizedGravityEPI {
public:
    PS::F64vec pos;
    PS::F64    eps2;
    void copyFromFP(const SPH & sph){ 
        pos  = sph.pos;
        eps2 = sph.ksr * sph.ksr * KernelSph::ksrhinv * KernelSph::ksrhinv;
    }
    void copyFromFP(const MassLess & msls){         
        pos  = msls.pos;
        eps2 = SPH::eps * SPH::eps;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
};

class SymmetrizedGravityEPJ {
public:
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    eps2;
    void copyFromFP(const SPH & sph){ 
        mass = sph.mass;
        pos  = sph.pos;
        eps2 = sph.ksr * sph.ksr * KernelSph::ksrhinv * KernelSph::ksrhinv;
    }
    void copyFromFP(const MassLess & msls){ 
        mass = msls.mass;
        pos  = msls.pos;
        eps2 = SPH::eps * SPH::eps;
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

template <class TGravityJ>
struct calcSymmetrizedGravityBasic {    

    void operator () (const SymmetrizedGravityEPI *epi,
                      const PS::S32 nip,
                      const TGravityJ *epj,
                      const PS::S32 njp,
                      Gravity *gravity) {

        for(PS::S32 i = 0; i < nip; i++) {
            PS::F64vec xi   = epi[i].pos;
            PS::F64    ei2  = epi[i].eps2;
            PS::F64vec ai   = 0.0;
            PS::F64    pi   = 0.0;
            PS::F64    etai = 0.0;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::F64vec xj  = epj[j].pos;
                PS::F64    ej2 = epj[j].eps2;
                PS::F64    mj  = epj[j].mass;
                
                PS::F64vec dx   = xi - xj;
                PS::F64    r2   = dx * dx;
                PS::F64    rei2 = r2 + ei2;
                PS::F64    rej2 = r2 + ej2;
                PS::F64    rei_inv = 1. / sqrt(rei2);
                PS::F64    rej_inv = 1. / sqrt(rej2);
                PS::F64    rei_inv3 = rei_inv * rei_inv * rei_inv;
                PS::F64    rej_inv3 = rej_inv * rej_inv * rej_inv;
                PS::F64vec dg2_i    = - 0.5 * mj * dx * (rei_inv3 + rej_inv3);
                PS::F64    dpi      = 0.5 * mj * (rei_inv + rej_inv);
                PS::F64    detai    = mj * rei_inv3;

                ai   += dg2_i;
                pi   -= dpi;
                etai += detai;
            }
            gravity[i].acc += ai;
            gravity[i].pot += pi;
            gravity[i].eta += etai;
        }
    }
};

template <class TGravityJ>
struct calcSymmetrizedGravityX86 {    

    void operator () (const SymmetrizedGravityEPI *epi,
                      const PS::S32 nip,
                      const TGravityJ *epj,
                      const PS::S32 njp,
                      Gravity *gravity) {

        v4df (*rsqrt)(v4df) = v4df::rsqrt_4th;
        PS::S32 nvector = v4df::getVectorLength();

        for(PS::S32 i = 0; i < nip; i += nvector) {
            v4df px_i(epi[i].pos[0], epi[i+1].pos[0], epi[i+2].pos[0], epi[i+3].pos[0]);
            v4df py_i(epi[i].pos[1], epi[i+1].pos[1], epi[i+2].pos[1], epi[i+3].pos[1]);
            v4df pz_i(epi[i].pos[2], epi[i+1].pos[2], epi[i+2].pos[2], epi[i+3].pos[2]);
            v4df e2_i(epi[i].eps2,   epi[i+1].eps2,   epi[i+2].eps2,   epi[i+3].eps2);
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

            PS::F64 buf_ax[nvector];
            PS::F64 buf_ay[nvector];
            PS::F64 buf_az[nvector];
            PS::F64 buf_pt[nvector];
            PS::F64 buf_et[nvector];

            ax_i.store(buf_ax);
            ay_i.store(buf_ay);
            az_i.store(buf_az);
            pt_i.store(buf_pt);
            et_i.store(buf_et);

            PS::S32 nii = ((nip - i) < nvector) ? (nip - i) : nvector;
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

#ifdef SYMMETRIZED_GRAVITY
#ifdef ENABLE_SIMDX86
typedef calcSymmetrizedGravityX86<SymmetrizedGravityEPJ>      calcGravityEPJ;
typedef calcSymmetrizedGravityX86<PS::SPJMonopoleSymmetrized> calcGravitySPJ;
#elif defined ENABLE_SIMDX86_SINGLE
typedef calcSymmetrizedGravitySingleX86<SymmetrizedGravityEPJ>      calcGravityEPJ;
typedef calcSymmetrizedGravitySingleX86<PS::SPJMonopoleSymmetrized> calcGravitySPJ;
#else
typedef calcSymmetrizedGravityBasic<SymmetrizedGravityEPJ>      calcGravityEPJ;
typedef calcSymmetrizedGravityBasic<PS::SPJMonopoleSymmetrized> calcGravitySPJ;
#endif
#endif
