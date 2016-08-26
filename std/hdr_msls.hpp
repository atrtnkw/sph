class MassLess : public SPH {
public:

    void copyFromForce(const Gravity & gravity);

    void writeAscii(FILE *fp) const {
        using namespace CodeUnit;
        PS::F64vec tpos = this->pos * UnitOfLength;
        PS::F64    tpot = this->pot * UnitOfEnergy;
        PS::F64vec tomg = RP::RotationalVelocity * UnitOfTimeInv;
        PS::F64vec tvec = tomg ^ tpos;
        fprintf(fp, " %+.16e %+.16e %+.16e", tpos[0], tpos[1], tpos[2]);
        fprintf(fp, " %+.16e %+.16e", tpot, tpot - 0.5 * (tvec * tvec));
        fprintf(fp, "\n");
    }

    void readHexa(FILE *fp) {
        PS::F64 (*cvt)(PS::U64) = convertU64ToF64;
        PS::U64 umass;
        PS::U64 upos[3];
        PS::U64 upot;
        fscanf(fp, "%llx", &umass);
        fscanf(fp, "%llx %llx %llx", &upos[0], &upos[1], &upos[2]);
        fscanf(fp, "%llx", &upot);
        this->mass   = cvt(umass);
        this->pos[0] = cvt(upos[0]);
        this->pos[1] = cvt(upos[1]);
        this->pos[2] = cvt(upos[2]);
        this->pot    = cvt(upot);
    }

    void writeHexa(FILE *fp) const {
        PS::U64 (*cvt)(PS::F64) = convertF64ToU64;
        fprintf(fp, "%llx", cvt(this->mass));
        fprintf(fp, " %llx %llx %llx", cvt(this->pos[0]), cvt(this->pos[1]), cvt(this->pos[2]));
        fprintf(fp, " %llx", cvt(this->pot));
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
    xmax = std::max(xmax, 0.); // ??? caution!!!
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
        PS::F64vec tv = RP::RotationalVelocity ^ msls[i].pos;
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

