namespace NthDimension {
#ifdef USE_IDEAL

    vndf (*calcVolumeInverse)(const vndf hi);
    PS::F64 (*calcPowerOfDimInverse)(const PS::F64 mass,
                                     const PS::F64 dens);
    PS::F64mat (*invertMatrix)(const PS::F64mat & tau);

    inline vndf calcVolumeInverse1D(const vndf hi) {
        return hi;
    }
    inline vndf calcVolumeInverse2D(const vndf hi) {
        return hi * hi;
    }
    inline vndf calcVolumeInverse3D(const vndf hi) {
        return hi * hi * hi;
    }

    inline PS::F64 calcPowerOfDimInverse1D(const PS::F64 mass,
                                           const PS::F64 dens) {
        return mass / dens;
    }
    inline PS::F64 calcPowerOfDimInverse2D(const PS::F64 mass,
                                           const PS::F64 dens) {
        return sqrt(mass / dens);
    }
    inline PS::F64 calcPowerOfDimInverse3D(const PS::F64 mass,
                                           const PS::F64 dens) {
        return pow(mass / dens, 1. / 3.);
    }

    inline PS::F64mat invertMatrix1D(const PS::F64mat & tau) {
        PS::F64mat c = 0.;
        c.xx = 1. / tau.xx;
        return c;
    }
    inline PS::F64mat invertMatrix2D(const PS::F64mat & tau) {
        PS::F64mat c = 0.;
        PS::F64 detinv = 1. / (tau.xx * tau.yy - tau.xy * tau.xy);
        c.xx =   tau.yy * detinv;
        c.yy =   tau.xx * detinv;
        c.xy = - tau.xy * detinv;
        return c;
    }
    inline PS::F64mat invertMatrix3D(const PS::F64mat & tau) {
        PS::F64mat c = 0.;
        PS::F64 detinv = 1. / (tau.xx * tau.yy * tau.zz
                               + tau.xy * tau.xz * tau.yz * 2.
                               - tau.xx * tau.yz * tau.yz
                               - tau.yy * tau.xz * tau.xz
                               - tau.zz * tau.xy * tau.xy);
        c.xx = (tau.yy * tau.zz - tau.yz * tau.yz) * detinv;
        c.yy = (tau.xx * tau.zz - tau.xz * tau.xz) * detinv;
        c.zz = (tau.xx * tau.yy - tau.xy * tau.xy) * detinv;
        c.xy = (tau.xz * tau.yz - tau.xy * tau.zz) * detinv;
        c.xz = (tau.xy * tau.yz - tau.xz * tau.yy) * detinv;
        c.yz = (tau.xz * tau.xy - tau.xx * tau.yz) * detinv;
        return c;
    }

    void setDimension(const PS::F64 ndim) {
        if(ndim == 1) {
            calcVolumeInverse     = calcVolumeInverse1D;
            calcPowerOfDimInverse = calcPowerOfDimInverse1D;
            invertMatrix          = invertMatrix1D;
            if(PS::Comm::getRank() == 0) {
                fprintf(stderr, "set 1D!\n");
            }
        } else if (ndim == 2) {        
            calcVolumeInverse     = calcVolumeInverse2D;
            calcPowerOfDimInverse = calcPowerOfDimInverse2D;
            invertMatrix          = invertMatrix2D;
            if(PS::Comm::getRank() == 0) {
                fprintf(stderr, "set 2D!\n");
            }
        } else {
            calcVolumeInverse     = calcVolumeInverse3D;
            calcPowerOfDimInverse = calcPowerOfDimInverse3D;
            invertMatrix          = invertMatrix3D;
            if(PS::Comm::getRank() == 0) {
                fprintf(stderr, "set 3D!\n");
            }
        }
    }
#elif defined USE_HELMHOLTZ

#ifdef FOR_TUBE_TEST

    inline vndf calcVolumeInverse(const vndf hi) {
        return hi;
    }

    inline PS::F64 calcPowerOfDimInverse(const PS::F64 mass,
                                         const PS::F64 dens) {
        return mass / dens;
    }

    inline PS::F64mat invertMatrix(const PS::F64mat & tau) {
        PS::F64mat c = 0.;
        c.xx = 1. / tau.xx;
        return c;
    }

    void setDimension(const PS::F64 ndim) {
        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "set 1D!\n");
        }
    }

#else


#ifdef USE_INTRINSICS

    inline vndf calcVolumeInverse(const vndf hi) {
        return hi * hi * hi;
    }

#else

    inline PS::F64 calcVolumeInverse(const PS::F64 hi) {
        return hi * hi * hi;
    }

#endif

    inline PS::F64 calcPowerOfDimInverse(const PS::F64 mass,
                                         const PS::F64 dens) {
        return pow(mass / dens, 1. / 3.);
    }

    inline PS::F64mat invertMatrix(const PS::F64mat & tau) {
        PS::F64mat c = 0.;
        PS::F64 detinv = 1. / (tau.xx * tau.yy * tau.zz
                               + tau.xy * tau.xz * tau.yz * 2.
                               - tau.xx * tau.yz * tau.yz
                               - tau.yy * tau.xz * tau.xz
                               - tau.zz * tau.xy * tau.xy);
        c.xx = (tau.yy * tau.zz - tau.yz * tau.yz) * detinv;
        c.yy = (tau.xx * tau.zz - tau.xz * tau.xz) * detinv;
        c.zz = (tau.xx * tau.yy - tau.xy * tau.xy) * detinv;
        c.xy = (tau.xz * tau.yz - tau.xy * tau.zz) * detinv;
        c.xz = (tau.xy * tau.yz - tau.xz * tau.yy) * detinv;
        c.yz = (tau.xz * tau.xy - tau.xx * tau.yz) * detinv;
        return c;
    }

    void setDimension(const PS::F64 ndim) {
        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "set 3D!\n");
        }
    }

#endif

#else
#error We have only two options: USE_IDEAL and USE_HELMHOLTZ.
#endif
    
}

namespace ND = NthDimension;
