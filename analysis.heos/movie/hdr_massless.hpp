class MassLess {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    ksr;
    PS::F64    dens;
    PS::F64    temp;
    PS::F64    abar;
    static PS::S64 nx;
    static PS::F64 xmin;
    static PS::F64 xmax;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {        
        this->pos = pos_new;
    }

    void copyFromForce(const Quantity & quantity) {
        this->dens = quantity.dens;
        this->temp = quantity.temp;
        this->abar = quantity.abar;
    }

    void writeAscii(FILE *fp) const {
        fprintf(fp, " %+e %+e", this->pos[0], this->pos[1]);
        fprintf(fp, " %+e %+e %+e", this->dens, this->temp, this->abar);
        fprintf(fp, "\n");
    }
};

PS::S64 MassLess::nx   = 512;
PS::F64 MassLess::xmin = -2e9;
PS::F64 MassLess::xmax = +2e9;

template <class Tmassless>
void generateMassLessParticle(Tmassless & msls) {
    PS::S64 nx    = MassLess::nx;
    PS::F64 xmin  = MassLess::xmin;
    PS::F64 xmax  = MassLess::xmax;
    PS::F64 dx    = (xmax - xmin) / (PS::F64)nx;
    PS::S64 nmsls = nx * nx;
    msls.setNumberOfParticleLocal(nmsls);
    for(PS::S64 i = 0; i < nmsls; i++) {
        msls[i].id     = -1;
        msls[i].mass   = 0.;
        msls[i].pos[0] = xmin + dx * (PS::F64)(i % nx);
        msls[i].pos[1] = xmin + dx * (PS::F64)(i / nx);
        msls[i].pos[2] = 0.;
        msls[i].ksr    = 4. * dx;
    }
}
