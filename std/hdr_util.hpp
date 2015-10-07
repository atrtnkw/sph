#pragma once

template <class Tsph>
void calcCenterOfMass(Tsph & sph,
                      PS::F64    & mc,
                      PS::F64vec & xc,
                      PS::F64vec & vc,
                      PS::S64 istar = -1) {
    PS::F64    mloc = 0.;
    PS::F64vec xloc = 0.;
    PS::F64vec vloc = 0.;

    assert(istar < 2);

    PS::S32 nloc = sph.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        if(sph[i].istar == istar || istar < 0) {
            mloc += sph[i].mass;
            xloc += sph[i].mass * sph[i].pos;
            vloc += sph[i].mass * sph[i].vel;
        }
    }

    mc = PS::Comm::getSum(mloc);
    xc = PS::Comm::getSum(xloc);
    vc = PS::Comm::getSum(vloc);
    PS::F64 minv = 1. / mc;
    xc *= minv;
    vc *= minv;    
}

template <class Tsph>
PS::F64 calcSystemSize(Tsph & sph) {
    PS::F64    m0, m1;
    PS::F64vec x0, x1, v0, v1;
    calcCenterOfMass(sph, m0, x0, v0, 0);
    calcCenterOfMass(sph, m1, x1, v1, 1);
    PS::F64vec dx = x0 - x1;
    PS::F64    dr = (m0 != 0. && m1 != 0.) ? sqrt(dx * dx)
        : std::numeric_limits<double>::max();
    return dr;
}

template <class Tsph>
void calcRotationalVelocity(Tsph & sph) {
    PS::F64    m0, m1;
    PS::F64vec x0, x1;
    PS::F64vec v0, v1;
    calcCenterOfMass(sph, m0, x0, v0, 0);
    calcCenterOfMass(sph, m1, x1, v1, 1);

    PS::F64vec axisv = x0 - x1;
    PS::F64    axis  = sqrt(axisv * axisv);
    PS::F64    vel   = sqrt(CodeUnit::GravityConstantInThisUnit * (m0 + m1) / axis);
    
    RP::RotationalVelocity[0] = 0.d;
    RP::RotationalVelocity[1] = 0.d;
    RP::RotationalVelocity[2] = vel / axis;
}

template <class Tsph,
          class Tmsls>
void reduceSeparation(Tsph & sph,
                      Tmsls & msls) {
//    static bool FirstStep    = true;
    static bool StopDamping2 = false;

//    if(FirstStep) {
//        FirstStep = false;
    if(RP::Time == 0.) {
        debugFunction("a");
        PS::F64    m1;
        PS::F64vec x1, v1;
        calcCenterOfMass(sph, m1, x1, v1, 1);

        PS::F64 r2maxloc = 0.d;
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].istar == 1) {
                PS::F64vec dr = sph[i].pos - x1;
                PS::F64    r2 = dr * dr;
                r2maxloc = std::max(r2, r2maxloc);
            }
        }        

        PS::F64 r2maxglb = PS::Comm::getMaxValue(r2maxloc);
        PS::F64 rmax     = sqrt(r2maxglb);
        PS::F64 rho2     = m1 / ((4.d * M_PI / 3.d * rmax * rmax * rmax));

        RP::ReductionTime    = 1. / (0.05 * sqrt(CodeUnit::GravityConstantInThisUnit * rho2));
        RP::ReductionTimeInv = 1. / RP::ReductionTime;
    }

    if(RP::Time - (PS::S64)(RP::Time / RP::MaximumTimestep) * RP::MaximumTimestep == 0.
       && RP::Time != 0.) {
        PS::F64    m0, m1;
        PS::F64vec x0, x1;
        PS::F64vec v0, v1;
        calcCenterOfMass(sph, m0, x0, v0, 0);
        calcCenterOfMass(sph, m1, x1, v1, 1);

        if(StopDamping2) {
            PS::F64 mc;
            PS::F64vec xc;
            PS::F64vec vc;
            calcCenterOfMass(sph, mc, xc, vc);            
            for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
                sph[i].pos -= xc;
                sph[i].vel -= vc;
                sph[i].vel += RP::RotationalVelocity ^ sph[i].pos;
            }    
            sph.writeParticleAscii("snap/final.dat");
            msls.writeParticleAscii("snap/msls_final.dat");
            PS::Finalize();
            exit(0);
        }

        PS::F64vec dx = x1 - x0;
        PS::F64    dr = sqrt(dx * dx);
        PS::F64    ds = dr / RP::ReductionTime * RP::MaximumTimestep;
        dx *= (- ds / dr);
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].istar == 1) {
                sph[i].pos += dx;
            }
        }
        PS::F64    mc;
        PS::F64vec xc;
        PS::F64vec vc;
        calcCenterOfMass(sph, mc, xc, vc);
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].pos -= xc;
        }

        PS::F64 xlag = searchLagrange1(msls, x0[0], x1[0]);
        bool StopDamping2Local = false;
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].istar == 0)
                continue;
            if(sph[i].pos[0] > xlag) {
                StopDamping2Local = true;
                break;
            }
        }
        StopDamping2 = PS::Comm::synchronizeConditionalBranchOR(StopDamping2Local);

        calcRotationalVelocity(sph);
   }
}

