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
        S32    id;
        F64    mass;
        F64vec pos;
        F64    eps2;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            id   = -1;
            mass = mom.getCharge();
            pos  = mom.getPos();
            eps2 = mom.getSofteningSquare();
        }
        void clear(){
            id   = 0;
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
    PS::S32    id;
    PS::F64vec pos;
    PS::F64    eps2;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        pos  = sph.pos;
        eps2 = sph.ksr * sph.ksr * KernelSph::ksrhinv * KernelSph::ksrhinv;
    }
    void copyFromFP(const MassLess & msls){         
        pos  = msls.pos;
        eps2 = 0.5 * SPH::eps * SPH::eps; // This is not accurate.
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
    PS::S32    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    eps2;
    void copyFromFP(const SPH & sph){ 
        id   = sph.id;
        mass = sph.mass;
        pos  = sph.pos;
        eps2 = sph.ksr * sph.ksr * KernelSph::ksrhinv * KernelSph::ksrhinv;
    }
    void copyFromFP(const MassLess & msls){ 
        mass = msls.mass;
        pos  = msls.pos;
        eps2 = 0.5 * SPH::eps * SPH::eps; // This is not accurate.
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
struct calcSymmetrizedGravity {    

    void operator () (const SymmetrizedGravityEPI *epi,
                      const PS::S32 nip,
                      const TGravityJ *epj,
                      const PS::S32 njp,
                      Gravity *gravity) {

        for(PS::S32 i = 0; i < nip; i++) {
            PS::S32    idi = epi[i].id;
            PS::F64vec xi  = epi[i].pos;
            PS::F64    ei2 = epi[i].eps2;
            PS::F64vec ai  = 0.0;
            PS::F64    pi  = 0.0;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S32    idj = epj[j].id;
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
                //PS::F64vec dg2_i    = - 0.5 * mj * dx * (rej_inv3 + rej_inv3);
                PS::F64vec dg2_i    = - 0.5 * mj * dx * (rei_inv3 + rej_inv3);
                PS::F64    dpi      = 0.5 * mj * (rei_inv + rej_inv);

                ai += dg2_i;
                pi -= dpi;
            }
            gravity[i].acc += ai;
            gravity[i].pot += pi;
        }
    }
};
