class MassLess {
public:
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    pot;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {        
        this->pos = pos_new;
    }

    void copyFromForce(const Gravity & gravity) {
        this->pot = CodeUnit::grav * (gravity.pot + this->mass / SPH::eps);
    }

    void writeAscii(FILE *fp) const {
        using namespace CodeUnit;
        PS::F64vec tpos  = this->pos  * UnitOfLength;
        PS::F64    tpot  = this->pot  * UnitOfEnergy;
        fprintf(fp, " %+.16e %+.16e %+.16e", tpos[0], tpos[1], tpos[2]);
        fprintf(fp, " %+.16e %+.16e", tpot, SPH::omg[2]);        
        fprintf(fp, "\n");
    }
};

template <class Tmassless,
          class Tsph>
void generateMassLessParticle(Tmassless & msls,
                              Tsph & sph) {
    PS::F64 xmaxloc = -1e60;
    PS::F64 xminloc = +1e60;
    PS::F64 hmaxloc = -1e60;
    PS::F64 hminloc = +1e60;
    PS::S32 nloc = sph.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        PS::F64 xi = sph[i].pos[0];
        PS::F64 hi = sph[i].ksr;
        xmaxloc = std::max(xi, xmaxloc);
        xminloc = std::min(xi, xminloc);
        hmaxloc = std::max(hi, hmaxloc);
        hminloc = std::min(hi, hminloc);
    }
    PS::F64 xmax = 2. * PS::Comm::getMaxValue(xmaxloc);
    PS::F64 xmin = 2. * PS::Comm::getMinValue(xminloc);
    PS::F64 hmax =      PS::Comm::getMaxValue(hmaxloc);
    PS::F64 hmin =      PS::Comm::getMinValue(hminloc);

    PS::F64 xwidth    = 0.25 * sqrt(hmax * hmin);
    PS::S32 nmassless = (PS::S32)((xmax - xmin) / xwidth);

    if(PS::Comm::getRank() == 0) {
        msls.setNumberOfParticleLocal(nmassless);
        for(PS::S32 i = 0; i < nmassless; i++) {
            MassLess ptcl;
            ptcl.mass   = 0.;
            ptcl.pos[0] = xmin + xwidth * (PS::F64)i;
            ptcl.pos[1] = 0.;
            ptcl.pos[2] = 0.;
            msls[i]     = ptcl;
        }
    }
}

template <class Tmassless>
PS::F64 searchLagrange1(Tmassless & msls,
                     PS::F64 x0,
                     PS::F64 x1) {
    MassLess ploc;
    PS::F64 potmaxloc = -1e60;
    PS::S32  nmsls = msls.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nmsls; i++) {
        PS::F64 xi = msls[i].pos[0];
        if(xi < x1 || x0 < xi)
            continue;
        PS::F64vec tv = SPH::omg ^ msls[i].pos;
        PS::F64 ipot = msls[i].pot - 0.5 * (tv * tv);
        if(ipot > potmaxloc) {
            ploc      = msls[i];
            potmaxloc = ipot;
        }
    }
    PS::S32 rloc = PS::Comm::getRank();
    PS::F64 potmaxglb;
    PS::S32 rglb;
    PS::Comm::getMaxValue(potmaxloc, rloc, potmaxglb, rglb);
    PS::F64 xlag = ploc.pos[0];
    PS::Comm::broadcast(&xlag, 1, rglb);
    return xlag;
}

