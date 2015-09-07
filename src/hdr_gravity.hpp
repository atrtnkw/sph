class GravityEPI {
public:
    PS::F64vec pos;
    void copyFromFP(const SPH & sph){ 
        pos  = sph.pos;
    }
    void copyFromFP(const MassLess & msls){         
        pos  = msls.pos;
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
    void copyFromFP(const SPH & sph){ 
        mass = sph.mass;
        pos  = sph.pos;
    }
    void copyFromFP(const MassLess & msls){ 
        mass = msls.mass;
        pos  = msls.pos;
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
};

namespace MyGlobal {
    const PS::S32 NumberOfIGroup = 64;
    const PS::S32 NumberOfJ      = 65536;
}

template <class TGravityJ>
struct calcGravity {    

    void operator () (const GravityEPI *epi,
                      const PS::S32 nip,
                      const TGravityJ *epj,
                      const PS::S32 njp,
                      Gravity *gravity) {

        using namespace MyGlobal;
        static __thread PS::F64 xi[NumberOfIGroup][3];
        static __thread PS::F64 ai[NumberOfIGroup][3];
        static __thread PS::F64 pi[NumberOfIGroup];
        static __thread PS::F64 xj[NumberOfJ][3];
        static __thread PS::F64 mj[NumberOfJ];

        assert(nip <= NumberOfIGroup);
        assert(njp <= NumberOfJ);

        for(PS::S32 i = 0; i < nip; i++) {
            xi[i][0] = epi[i].pos[0];
            xi[i][1] = epi[i].pos[1];
            xi[i][2] = epi[i].pos[2];
            ai[i][0] = 0.0;
            ai[i][1] = 0.0;
            ai[i][2] = 0.0;
            pi[i]    = 0.0;
        }

        for(PS::S32 j = 0; j < njp; j++) {
            xj[j][0] = epj[j].pos[0];
            xj[j][1] = epj[j].pos[1];
            xj[j][2] = epj[j].pos[2];
            mj[j]    = epj[j].mass;
        }

        PS::S32 devid = PS::Comm::getThreadNum();

        g5_set_xmjMC(devid, 0, njp, xj, mj);
        g5_set_nMC(devid, njp);
        g5_calculate_force_on_xMC(devid, xi, ai, pi, nip);

        for(PS::S32 i = 0; i < nip; i++) {
            gravity[i].acc[0] += ai[i][0];
            gravity[i].acc[1] += ai[i][1];
            gravity[i].acc[2] += ai[i][2];
            gravity[i].pot    -= pi[i];
        }        
    }
};

#if not defined SYMMETRIZED_GRAVITY
typedef calcGravity<GravityEPJ>      calcGravityEPJ;
typedef calcGravity<PS::SPJMonopole> calcGravitySPJ;
#endif
