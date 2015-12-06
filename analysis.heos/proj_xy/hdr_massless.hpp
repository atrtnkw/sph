class MassLess {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    ksr;
    PS::F64    dens;
    PS::F64    temp;
    PS::F64    abar;
    PS::F64    shck;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {        
        this->pos = pos_new;
    }

    void copyFromForce(const Quantity & quantity);

    void writeAscii(FILE *fp) const {
        using namespace CodeUnit;
        PS::F64vec tpos  = this->pos  * UnitOfLength;
        PS::F64    tdens = this->dens * UnitOfDensity;
        fprintf(fp, " %+e %+e", tpos[0], tpos[1]);
        fprintf(fp, " %+e %+e %+e %+e", tdens, this->temp, this->abar, this->shck);
        fprintf(fp, "\n");
    }
};

template <class Tmassless>
void generateMassLessParticle(Header & hdr,
                              Tmassless & msls) {
    PS::S64 nx    = hdr.msls_nx;
    PS::F64 xmin  = hdr.msls_xmin;
    PS::F64 xmax  = hdr.msls_xmax;
    PS::F64 dx    = (xmax - xmin) / (PS::F64)nx;
    PS::S64 nmsls = nx * nx;
    msls.setNumberOfParticleLocal(nmsls);
    for(PS::S64 i = 0; i < nmsls; i++) {
        msls[i].id     = -1;
        msls[i].mass   = 0.;
        msls[i].pos[0] = xmin + dx * (PS::F64)(i % nx);
        msls[i].pos[1] = xmin + dx * (PS::F64)(i / nx);
        //msls[i].pos[2] = 0.;
        msls[i].pos[2] = hdr.cpos[2];
        msls[i].ksr    = 4. * dx;
    }
}
