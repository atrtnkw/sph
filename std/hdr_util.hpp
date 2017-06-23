#pragma once

inline PS::F64 getRandomNumber() {
    return ((PS::F64) rand() / ((PS::F64)RAND_MAX + 1.));
}

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
void calcCenterOfMassOfAcceleration(Tsph & sph,
                                    PS::F64    & mc,
                                    PS::F64vec & ac,
                                    PS::S64 istar = -1) {
    PS::F64    mloc = 0.;
    PS::F64vec aloc = 0.;

    assert(istar < 2);

    PS::S32 nloc = sph.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        if(sph[i].istar == istar || istar < 0) {
            mloc += sph[i].mass;
            aloc += sph[i].mass * sph[i].acc;
        }
    }

    mc = PS::Comm::getSum(mloc);
    ac = PS::Comm::getSum(aloc);
    PS::F64 minv = 1. / mc;
    ac *= minv;
}

template <class Tsph,
          class Tbhns>
void calcCenterOfMass(Tsph & sph,
                      Tbhns & bhns,
                      PS::F64    & mc,
                      PS::F64vec & xc,
                      PS::F64vec & vc) {
    PS::F64 m0, m1;
    PS::F64vec x0, x1;
    PS::F64vec v0, v1;                
    calcCenterOfMass(bhns, m0, x0, v0);
    calcCenterOfMass(sph,  m1, x1, v1);
    mc = m0 + m1;
    xc = (m0 * x0 + m1 * x1) / (m0 + m1);
    vc = (m0 * v0 + m1 * v1) / (m0 + m1);
}

template <class Tsph,
          class Tbhns>
void calcCenterOfMassOfAcceleration(Tsph & sph,
                                    Tbhns & bhns,
                                    PS::F64    & mc,
                                    PS::F64vec & xc,
                                    PS::F64vec & vc,
                                    PS::F64vec & ac) {
    calcCenterOfMass(sph, bhns, mc, xc, vc);
    PS::F64    m0, m1;
    PS::F64vec a0, a1;
    calcCenterOfMassOfAcceleration(bhns, m0, a0);
    calcCenterOfMassOfAcceleration(sph,  m1, a1);
    ac = (m0 * a0 + m1 * a1) / (m0 + m1);
}

template <class Tsph,
          class Tbhns>
PS::F64 calcSystemSize(Tsph & sph,
                       Tbhns & bhns) {
    PS::F64    m0, m1;
    PS::F64vec x0, x1, v0, v1;
    if(RP::FlagBinary == 0) {
        calcCenterOfMass(sph, m0, x0, v0, 0);
    } else if(RP::FlagBinary == 1) {
        calcCenterOfMass(bhns, m0, x0, v0, 0);
    } else if(RP::FlagBinary == 2) {
        calcCenterOfMass(sph, m0, x0, v0, 0);
    } else {
        ;
    }
    calcCenterOfMass(sph, m1, x1, v1, 1);
    PS::F64vec dx = x0 - x1;
    PS::F64    dr = (m0 != 0. && m1 != 0.) ? sqrt(dx * dx)
        : std::numeric_limits<double>::max();
    return dr;
}

template <class Tsph,
          class Tbhns>
void calcRotationalVelocity(Tsph & sph,
                            Tbhns & bhns) {
    PS::F64    m0, m1;
    PS::F64vec x0, x1;
    PS::F64vec v0, v1;
    if(RP::FlagBinary == 0) {
        calcCenterOfMass(sph, m0, x0, v0, 0);
    } else if(RP::FlagBinary == 1) {
        calcCenterOfMass(bhns, m0, x0, v0, 0);
    } else {
        ;
    }
    calcCenterOfMass(sph, m1, x1, v1, 1);

    PS::F64vec axisv = x0 - x1;
    PS::F64    axis  = sqrt(axisv * axisv);
    PS::F64    vel   = sqrt(CodeUnit::GravityConstantInThisUnit * (m0 + m1) / axis);
    
    RP::RotationalVelocity[0] = 0.;
    RP::RotationalVelocity[1] = 0.;
    RP::RotationalVelocity[2] = vel / axis;
}

template <class Tsph,
          class Tbhns,
          class Tmsls>
void reduceSeparation(Tsph & sph,
                      Tbhns & bhns,
                      Tmsls & msls) {
    static bool StopDamping2 = false;

    if(RP::Time - (PS::S64)(RP::Time / RP::MaximumTimestep) * RP::MaximumTimestep == 0.) {
        PS::F64    m1;
        PS::F64vec x1, v1;
        calcCenterOfMass(sph, m1, x1, v1, 1);

        PS::F64 r2maxloc = 0.;
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].istar == 1) {
                PS::F64vec dr = sph[i].pos - x1;
                PS::F64    r2 = dr * dr;
                r2maxloc = std::max(r2, r2maxloc);
            }
        }        

        PS::F64 r2maxglb = PS::Comm::getMaxValue(r2maxloc);
        PS::F64 rmax     = sqrt(r2maxglb);
        PS::F64 rho2     = m1 / ((4. * M_PI / 3. * rmax * rmax * rmax));

        RP::ReductionTime    = 1. / (0.05 * sqrt(CodeUnit::GravityConstantInThisUnit * rho2));
        RP::ReductionTimeInv = 1. / RP::ReductionTime;
    }

    if(RP::Time - (PS::S64)(RP::Time / RP::MaximumTimestep) * RP::MaximumTimestep == 0.
       && RP::Time != 0.) {
        PS::F64    m0, m1;
        PS::F64vec x0, x1;
        PS::F64vec v0, v1;
        if(RP::FlagBinary == 0) {
            calcCenterOfMass(sph, m0, x0, v0, 0);
        } else if (RP::FlagBinary == 1) {
            calcCenterOfMass(bhns, m0, x0, v0, 0);
        } else {
            ;
        }
        calcCenterOfMass(sph, m1, x1, v1, 1);

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

        if(StopDamping2) {
            PS::F64 mc;
            PS::F64vec xc;
            PS::F64vec vc;
            if(RP::FlagBinary == 0) {
                calcCenterOfMass(sph, mc, xc, vc);
            } else if(RP::FlagBinary == 1) {
                calcCenterOfMass(sph, bhns, mc, xc, vc);
            } else {
                ;
            }
            for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
                sph[i].pos -= xc;
                sph[i].vel -= vc;
                sph[i].vel += RP::RotationalVelocity ^ sph[i].pos;
            }    
            for(PS::S32 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
                bhns[i].pos -= xc;
                bhns[i].vel -= vc;
                bhns[i].vel += RP::RotationalVelocity ^ bhns[i].pos;
            }    
            if(RP::FlagDivideFile == 0) {
                sph.writeParticleAscii("snap/final.dat");
            } else {
                char filename[64];
                sprintf(filename, "snap/final");
                sph.writeParticleAscii(filename, "%s_p%06d_i%06d.dat");
            }
            msls.writeParticleAscii("snap/msls_final.dat");
            bhns.writeParticleAscii("snap/bhns_final.dat");
            PS::Finalize();
            exit(0);
        }

        PS::F64vec dx = x1 - x0;
        PS::F64    dr = sqrt(dx * dx);
        PS::F64    ds = dr / RP::ReductionTime * RP::MaximumTimestep; // ds: Delta r
        dx *= (- ds / dr);  // dx: Delta x
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].istar == 1) {
                sph[i].pos += dx;
            }
        }
        PS::F64    mc;
        PS::F64vec xc;
        PS::F64vec vc;
        if(RP::FlagBinary == 0) {
            calcCenterOfMass(sph, mc, xc, vc);
        } else if (RP::FlagBinary == 1) {
            calcCenterOfMass(sph, bhns, mc, xc, vc);
        } else {
            ;
        }
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].pos -= xc;
        }
        for(PS::S32 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
            bhns[i].pos -= xc;
        }
        calcRotationalVelocity(sph, bhns);
        {
            PS::F64    m0, m1;
            PS::F64vec x0, x1;
            PS::F64vec v0, v1;
            if(RP::FlagBinary == 0) {
                calcCenterOfMass(sph, m0, x0, v0, 0);
            } else if (RP::FlagBinary == 1) {
                calcCenterOfMass(bhns, m0, x0, v0, 0);
            } else {
                ;
            }
            calcCenterOfMass(sph, m1, x1, v1, 1);
            RP::CenterOfMassOfStar0 = x0;
            RP::CenterOfMassOfStar1 = x1;
        }

    }
}

template <class Tbhns>
void broadcastBlackHoleNeutronStar(const Tbhns & bhns,
                                   PS::F64    & mass,
                                   PS::F64vec & pos,
                                   PS::F64vec & vel,
                                   PS::F64    & size) {
    PS::S32    nbhns = bhns.getNumberOfParticleLocal();
    PS::F64    mloc  = 0.;
    PS::F64vec ploc  = 0.;
    PS::F64vec vloc  = 0.;
    PS::F64    sloc  = 0.;
    if(nbhns == 1) {
        mloc = bhns[0].mass;
        ploc = bhns[0].pos;
        vloc = bhns[0].vel;
        sloc = bhns[0].size;
    } else if(nbhns == 0) {
        mloc = 0.;
        ploc = 0.;
        vloc = 0.;
        sloc = 0.;
    } else {
        fprintf(stderr, "Not supported in this case (in function %s)!\n", __FUNCTION__);
        PS::Abort();
    }
    mass = PS::Comm::getSum(mloc);
    pos  = PS::Comm::getSum(ploc);
    vel  = PS::Comm::getSum(vloc);
    size = PS::Comm::getSum(sloc);
}

template <class Tbhns>
void correctBlackHoleNeutronStar(Tbhns & bhns,
                                 const PS::F64    mass,
                                 const PS::F64vec masspos,
                                 const PS::F64vec momentum) {
    PS::S32    nbhns = bhns.getNumberOfParticleLocal();
    if(nbhns == 1) {
        bhns[0].pos  = (bhns[0].mass * bhns[0].pos + masspos)  / (bhns[0].mass + mass);
        bhns[0].vel  = (bhns[0].mass * bhns[0].vel + momentum) / (bhns[0].mass + mass);
        bhns[0].mass = bhns[0].mass + mass;
    } else if(nbhns == 0) {
        ;
    } else {
        fprintf(stderr, "Not supported in this case (in function %s)!\n", __FUNCTION__);
        PS::Abort();
    }    
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls,
          class Tdensity,
          class Thydro,
          class Tgravity>
void absorbParticleIntoBlackHoleNeutronStar(Tdinfo & dinfo,
                                            Tsph & sph,
                                            Tbhns & bhns,
                                            Tmsls & msls,
                                            Tdensity & density,
                                            Thydro & hydro,
                                            Tgravity & gravity) {
    PS::F64    mbhns;
    PS::F64vec pbhns, vbhns;
    PS::F64    sbhns;
    broadcastBlackHoleNeutronStar(bhns, mbhns, pbhns, vbhns, sbhns);

    PS::F64    ekloc = 0.;
    PS::F64    ephi0 = calcPotentialEnergy(sph, bhns);
    PS::F64    uloc0 = 0.;
    PS::F64    gloc0 = 0.;
    PS::S32    ndloc = 0;
    PS::F64    msloc = 0.;
    PS::F64vec mxloc = 0.;
    PS::F64vec mmloc = 0.;
    PS::F64    etot0 = calcEnergy(sph, bhns);

    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); ) {
        PS::F64vec dx = sph[i].pos - pbhns;
        PS::F64    r2 = dx * dx;
        PS::F64vec dv = sph[i].vel - vbhns;
        // tanikawa modified this 170623 FROM
        //if(r2 >= sph[i].ksr * sph[i].ksr) {
        if(r2 >= (sph[i].ksr + sbhns) * (sph[i].ksr + sbhns)) {
        // tanikawa modified this 170623 TO
            i++;
            continue;
        }        
        fprintf(RP::FilePointerForDebug, "\n");
        fprintf(RP::FilePointerForDebug, "### Absorbed %16.10f\n", RP::Time * CodeUnit::UnitOfTime);
        fprintf(RP::FilePointerForDebug, "### ");
        sph[i].writeAscii(RP::FilePointerForDebug);
        fprintf(RP::FilePointerForDebug, "### dx: %+e %+e %+e", dx[0] * CodeUnit::UnitOfLength,
                dx[1] * CodeUnit::UnitOfLength, dx[2] * CodeUnit::UnitOfLength);
        fprintf(RP::FilePointerForDebug, " dv: %+e %+e %+e\n", dv[0] * CodeUnit::UnitOfVelocity,
                dv[1] * CodeUnit::UnitOfVelocity, dv[2] * CodeUnit::UnitOfVelocity);
        fflush(RP::FilePointerForDebug);
        ndloc++;
        ekloc += 0.5 * sph[i].mass * (sph[i].vel * sph[i].vel);
        uloc0 +=       sph[i].mass *  sph[i].uene;
        gloc0 +=       sph[i].mass *  sph[i].gcor3;
        msloc +=       sph[i].mass;
        mxloc +=       sph[i].mass *  sph[i].pos;
        mmloc +=       sph[i].mass *  sph[i].vel;
        PS::S32 nsph = sph.getNumberOfParticleLocal();
        sph[i] = sph[nsph-1];
        sph.setNumberOfParticleLocal(nsph-1);
    }
    PS::S32    ndglb = PS::Comm::getSum(ndloc);
    PS::F64    kglb0 = PS::Comm::getSum(ekloc) + 0.5 * mbhns * (vbhns * vbhns);
    PS::F64    uglb0 = PS::Comm::getSum(uloc0);
    PS::F64    gglb0 = PS::Comm::getSum(gloc0);
    PS::F64    msglb = PS::Comm::getSum(msloc);
    PS::F64vec mxglb = PS::Comm::getSum(mxloc);
    PS::F64vec mmglb = PS::Comm::getSum(mmloc);

    if(ndglb > 0) {
        if(RP::FlagBinary == 1) {
            correctBlackHoleNeutronStar(bhns, msglb, mxglb, mmglb);
        } else if (RP::FlagBinary == 2) {
            // tanikawa add this 170623 FROM
            bhns.setNumberOfParticleLocal(0);
            // tanikawa add this 170623 TO
        } else {
            fprintf(stderr, "Something wrong happens!\n");
            PS::Abort();
        }

        calcSPHKernel(dinfo, sph, bhns, msls, density, hydro, gravity);

        PS::F64 etot1 = calcEnergy(sph, bhns);

        broadcastBlackHoleNeutronStar(bhns, mbhns, pbhns, vbhns, sbhns);
        PS::F64 kglb1 = 0.5 * mbhns * (vbhns * vbhns);
        PS::F64 ephi1 = calcPotentialEnergy(sph, bhns);
        RP::AbsorbedEnergyTotal += (kglb1 + ephi1) - (kglb0 + uglb0 + ephi0 - gglb0);

    }
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls,
          class Tdensity,
          class Thydro,
          class Tgravity>
void absorbParticleIntoIntermediateMassBlackHole(Tdinfo & dinfo,
                                                 Tsph & sph,
                                                 Tbhns & bhns,
                                                 Tmsls & msls,
                                                 Tdensity & density,
                                                 Thydro & hydro,
                                                 Tgravity & gravity) {
    static __thread Tbhns imbh;
    static __thread bool first = true;
    if(first) {
        imbh.initialize();
        imbh.createParticle(0);
        first = false;
    }
    if(PS::Comm::getRank() == 0) {
        using namespace CodeUnit;
        imbh.setNumberOfParticleLocal(1);
        imbh[0].mass = 0.;
        imbh[0].pos  = 0.;
        imbh[0].vel  = 0.;
        imbh[0].pot  = 0.;
        imbh[0].eps  = 1.;
        imbh[0].size = GravityConstantInThisUnit * BlackHoleMassInThisUnit /
            (SpeedOfLightInThisUnit * SpeedOfLightInThisUnit);
    } else {
        imbh.setNumberOfParticleLocal(0);
    }
    absorbParticleIntoBlackHoleNeutronStar(dinfo, sph, imbh, msls, density, hydro, gravity);
}

template <class Tsph>
PS::F64 calcDetonationVelocity(Tsph & sph) {
    static PS::F64 vdet = 0.;
#ifdef FOR_TUBE_TEST
    static PS::F64 xdet_prev = 0.;
    static PS::F64 time_prev = 0.;
    static PS::F64 tout      = 0.;
    static PS::F64 dtout     = RP::TimestepAscii / 4.;
    if(RP::Time >= tout) {
        PS::F64 tmax_loc = 0.;
        PS::F64 xdet_loc = 0.;
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].pos[0] < 0.) {
                continue;
            }
            if(sph[i].temp > tmax_loc) {
                tmax_loc = sph[i].temp;
                xdet_loc = sph[i].pos[0];
            }
        }
        PS::F64 tmax_glb;
        PS::S32 irank;
        PS::Comm::getMaxValue(tmax_loc, PS::Comm::getRank(), tmax_glb, irank);
        PS::F64 xdet_glb = xdet_loc;
        PS::Comm::broadcast(&xdet_glb, 1, irank);
        vdet = ((xdet_glb - xdet_prev) / (RP::Time - time_prev)) * CodeUnit::UnitOfVelocity;
        xdet_prev = xdet_glb;
        time_prev = RP::Time;
        tout     += dtout;
    }
#else
    //assert(NULL);
#endif
    return vdet;
}

template <class Tsph>
PS::F64 calcKernelSupportRadiusMaximum(Tsph & sph) {
    PS::F64 ksrmax = RP::KernelSupportRadiusMaximum;
    if(RP::FlagBinary == 0) {
//        if(RP::FlagDamping == 4) {
        PS::F64 mc;
        PS::F64vec xc, vc;
        calcCenterOfMass(sph, mc, xc, vc, 0);
        PS::F64 r2sumloc = 0.;
        PS::F64 dnsumloc = 0.;
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            PS::F64vec dx = sph[i].pos - xc;
            r2sumloc += (sph[i].dens * (dx * dx));
            dnsumloc +=  sph[i].dens;
        }
        PS::F64 r2sumglb = PS::Comm::getSum(r2sumloc);
        PS::F64 dnsumglb = PS::Comm::getSum(dnsumloc);
        ksrmax = sqrt(r2sumglb / dnsumglb);
//        }
    } else if(RP::FlagBinary == 1) {
        PS::F64 mc;
        PS::F64vec xc, vc;
        calcCenterOfMass(sph, mc, xc, vc, 1);
#if 0
        PS::F64 r2sumloc = 0.;
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            PS::F64vec dx = sph[i].pos - xc;
            r2sumloc += (dx * dx);
        }
        PS::F64 r2sumglb = PS::Comm::getSum(r2sumloc);
        ksrmax = sqrt(r2sumglb / (PS::F64)sph.getNumberOfParticleGlobal());
#else
        PS::F64 r2sumloc = 0.;
        PS::F64 dnsumloc = 0.;
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            PS::F64vec dx = sph[i].pos - xc;
            r2sumloc += (sph[i].dens * (dx * dx));
            dnsumloc +=  sph[i].dens;
        }
        PS::F64 r2sumglb = PS::Comm::getSum(r2sumloc);
        PS::F64 dnsumglb = PS::Comm::getSum(dnsumloc);
        ksrmax = sqrt(r2sumglb / dnsumglb);
#endif
    } else if(RP::FlagBinary == 2) {
        ;
    } else {
        assert(NULL);
    }
    return ksrmax;
}

template <class TParticle>
inline PS::F64mat calcQuadrupoleMomentDot2Component(const TParticle & ptcl,
                                                    const PS::F64vec & xc,
                                                    const PS::F64vec & vc,
                                                    const PS::F64vec & ac) {
    PS::F64mat qq = 0.;
    PS::F64vec dx = ptcl.pos - xc;
    PS::F64vec dv = ptcl.vel - vc;
    PS::F64vec da = ptcl.acc - ac;
    qq.xx = ptcl.mass * (2. * dv[0] * dv[0] + dx[0] * da[0] + da[0] * dx[0]);
    qq.xy = ptcl.mass * (2. * dv[0] * dv[1] + dx[0] * da[1] + da[0] * dx[1]);
    qq.xz = ptcl.mass * (2. * dv[0] * dv[2] + dx[0] * da[2] + da[0] * dx[2]);
    qq.yy = ptcl.mass * (2. * dv[1] * dv[1] + dx[1] * da[1] + da[1] * dx[1]);
    qq.yz = ptcl.mass * (2. * dv[1] * dv[2] + dx[1] * da[2] + da[1] * dx[2]);
    qq.zz = ptcl.mass * (2. * dv[2] * dv[2] + dx[2] * da[2] + da[2] * dx[2]);
    return qq;
}

template <class Tsph,
          class Tbhns>
void calcQuadrupoleMomentDot2(Tsph & sph,
                              Tbhns & bhns) {
    if(RP::FlagDamping != 0) {
        return;
    }

    if(RP::FlagBinary == 0) {
        ;
    } else if(RP::FlagBinary == 1) {
        PS::F64    mc;
        PS::F64vec xc, vc, ac;
        calcCenterOfMassOfAcceleration(sph, bhns, mc, xc, vc, ac);
        PS::F64mat qloc = 0.;
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            qloc = qloc + calcQuadrupoleMomentDot2Component(sph[i], xc, vc, ac);
        }
        for(PS::S64 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
            qloc = qloc + calcQuadrupoleMomentDot2Component(bhns[i], xc, vc, ac);
        }
        PS::F64mat qglb = 0.;
        qglb.xx = PS::Comm::getSum(qloc.xx);
        qglb.xy = PS::Comm::getSum(qloc.xy);
        qglb.xz = PS::Comm::getSum(qloc.xz);
        qglb.yy = PS::Comm::getSum(qloc.yy);
        qglb.yz = PS::Comm::getSum(qloc.yz);
        qglb.zz = PS::Comm::getSum(qloc.zz);
        if(PS::Comm::getRank() == 0) {
            using namespace CodeUnit;
            PS::F64 cc = (UnitOfMass * UnitOfLength * UnitOfLength * UnitOfTimeInv * UnitOfTimeInv)
                * (GravityConstant / (SpeedOfLight * SpeedOfLight * SpeedOfLight * SpeedOfLight));
            qglb.xx *= cc;
            qglb.xy *= cc;
            qglb.xz *= cc;
            qglb.yy *= cc;
            qglb.yz *= cc;
            qglb.zz *= cc;
            FILE *fp = RP::FilePointerForQuad;
            fprintf(fp, "%16.10f %+.16e %+.16e %+.16e %+.16e %+.16e %+.16e\n", RP::Time,
                    qglb.xx, qglb.xy, qglb.xz, qglb.yy, qglb.yz, qglb.zz);
            fflush(fp);
        }
    } else if(RP::FlagBinary == 2) {
        ;
    } else {
        assert(NULL);
    }
}
