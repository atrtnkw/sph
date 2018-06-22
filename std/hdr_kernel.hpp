#pragma once

namespace SmoothingKernel {
#ifdef USE_IDEAL
    PS::F64 eta;
    PS::F64 ksrh;
    PS::F64 ksrhinv;
    PS::F64 dim;
    PS::F64 ceff0;
    PS::F64 ceff1;
    PS::F64 ceff0x;

    vndf (*kernel0th)(const vndf q);
    vndf (*kernel1st)(const vndf q);

    inline vndf calcCubicSpline0th(const vndf q) {
        vndf t1 = vndf(1.)  - q;
        vndf t2 = vndf(0.5) - q;
        t1 = vndf::max(t1, 0.);
        t2 = vndf::max(t2, 0.);
        vndf t13 = t1 * t1 * t1;
        vndf t23 = t2 * t2 * t2;
        vndf w0  = t13 - vndf(4.) * t23;
        w0 *= vndf(ceff0x);
        
        return w0;
    }
    inline vndf calcCubicSpline1st(const vndf q) {
        vndf t1 = vndf(1.)  - q;
        vndf t2 = vndf(0.5) - q;
        vndf t3 = vndf(0.33333333333333333) - q;
        t1 = vndf::max(t1, 0.);
        t2 = vndf::max(t2, 0.);
        t3 = vndf::max(t3, 0.);
        vndf t12 = t1 * t1;
        vndf t22 = t2 * t2;
        vndf t32 = t3 * t3;
        vndf w0  = vndf::nmadd(t12, vndf(4.), t22);
        w0  = vndf::madd(w0, vndf(3.), t32);
        w0 *= vndf(ceff1);
        return w0;
    }

    inline vndf calcWendlandC20th1D(const vndf r) {
        vndf rmin = vndf::max(vndf(1.) - r, 0.);
        return vndf(ceff0) * rmin * rmin * rmin * (vndf(1.) + vndf(3.) * r);
    }
    inline vndf calcWendlandC21st1D(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin * rmin;
        vndf rmin3 = rmin * rmin2;
        return vndf(ceff1) * (rmin3 - rmin2 * (vndf(1.) + vndf(3.) * r));
    }

    inline vndf calcWendlandC20thMD(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin * rmin;
        return vndf(ceff0) * rmin2 * rmin2 * (vndf(1.) + vndf(4.) * r);
    }
    inline vndf calcWendlandC21stMD(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin  * rmin;
        vndf rmin3 = rmin  * rmin2;
        vndf rmin4 = rmin2 * rmin2;
        return vndf(ceff1) * (rmin4 - rmin3 * (vndf(1.) + vndf(4.) * r));
    }

    void setKernel(const KernelType kn,
                   const PS::F64 ndim) {
        if(kn == CubicSpline) {
            eta       = 1.2;
            ksrh      = 2.0;
            ksrhinv   = 1. / ksrh;
            kernel0th = calcCubicSpline0th;
            kernel1st = calcCubicSpline1st;
            if(ndim == 1) {
                dim   = +1.0;
                ceff0 = +1.3333333333333333;
                ceff1 = -8.;
                if(PS::Comm::getRank() == 0) {
                    fprintf(stderr, "set 1D cubic spline!\n");
                }
            } else if (ndim == 2) {
                dim   = +2.0;
                ceff0 = +1.8189136353359468;
                ceff1 = -10.913481812015681;
                if(PS::Comm::getRank() == 0) {
                    fprintf(stderr, "set 2D cubic spline!\n");
                }
            } else {
                dim   = +3.0;
                ceff0  = +2.5464790894703255;
                ceff1  = -15.278874536821953;
                if(PS::Comm::getRank() == 0) {
                    fprintf(stderr, "set 3D cubic spline!\n");
                }
            }
            ceff0x = ceff0 * 2.;
        } else if(kn == WendlandC2) {
            eta = 1.6;
            dim = (PS::F64)ndim;
            if(ndim == 1) {
                ksrh    = 1.620185;
                ceff0   = +1.25;
                ceff1   = +3.75;
                kernel0th = calcWendlandC20th1D;
                kernel1st = calcWendlandC21st1D;
                if(PS::Comm::getRank() == 0) {
                    fprintf(stderr, "set 1D WendlandC2!\n");
                }
            } else if(ndim == 2) {
                ksrh    = 1.897367;
                ceff0   = +2.228169203286535005e+00;
                ceff1   = +8.912676813146140020e+00;
                kernel0th = calcWendlandC20thMD;
                kernel1st = calcWendlandC21stMD;
                if(PS::Comm::getRank() == 0) {
                    fprintf(stderr, "set 2D WendlandC2!\n");
                }
            } else {
                ksrh    = 1.936492;
                ceff0   = +3.342253804929802286e+00;
                ceff1   = +1.336901521971920914e+01;
                kernel0th = calcWendlandC20thMD;
                kernel1st = calcWendlandC21stMD;
                if(PS::Comm::getRank() == 0) {
                    fprintf(stderr, "set 3D WendlandC2!\n");
                }
            }
            ksrhinv = 1. / ksrh;
        }
    }
#elif defined USE_HELMHOLTZ

#ifdef FOR_TUBE_TEST

#if 1
    const PS::F64 eta     = 1.6;
    //const PS::F64 eta     = 3.2;
    const PS::F64 ksrh    = 1.620185;
    const PS::F64 ksrhinv = 1. / ksrh;
    const PS::F64 dim     = 1.;
    const PS::F64 ceff0   = +1.25;
    const PS::F64 ceff1   = +3.75;

    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = ((1. - r > 0.) ? 1. - r : 0.);
        return ceff0 * rmin * rmin * rmin * (1. + 3. * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = ((1. - r > 0.) ? 1. - r : 0.);
        PS::F64 rmin2 = rmin * rmin;
        PS::F64 rmin3 = rmin * rmin2;
        return ceff1 * (rmin3 - rmin2 * (1. + 3. * r));
    }

    inline vndf kernel0th(const vndf r) {
        vndf rmin = vndf::max(vndf(1.) - r, 0.);
        return vndf(ceff0) * rmin * rmin * rmin * (vndf(1.) + vndf(3.) * r);
    }
    inline vndf kernel1st(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin * rmin;
        vndf rmin3 = rmin * rmin2;
        return vndf(ceff1) * (rmin3 - rmin2 * (vndf(1.) + vndf(3.) * r));
    }
    void setKernel(const KernelType kn,
                   const PS::F64 ndim) {
        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "set 1D WendlandC2!\n");
        }
    };
#else
    const  PS::F64 eta   = 1.6;
    const  PS::F64 dim   = +1.0;
    const  PS::F64 ksrh  = +1.936492;
    const  PS::F64 ksrhinv = 1. / ksrh;
    const  PS::F64 ceff0 = +1.5;

    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = ((1. - r > 0.) ? 1. - r : 0.);
        PS::F64 rmin2 = rmin * rmin;
        return ceff0 * rmin2 * rmin2 * rmin * (1. + 5. * r + 8. * r * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = ((1. - r > 0.) ? 1. - r : 0.);
        PS::F64 rmin2 = rmin  * rmin;
        PS::F64 rmin4 = rmin2 * rmin2;
        return ceff0 * rmin4 * (rmin * (5. + 16. * r) - (5. + 25. * r + 40 * r * r));
    }

    inline vndf kernel0th(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin * rmin;
        return vndf(ceff0) * rmin2 * rmin2 * rmin
            * (vndf(1.) + vndf(5.) * r + vndf(8.) * r * r);
    }
    inline vndf kernel1st(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin  * rmin;
        vndf rmin4 = rmin2 * rmin2;
        return vndf(ceff0) * rmin4 * (rmin * (vndf(5.) + vndf(16.) * r)
                                      - (vndf(5.) + vndf(25.) * r + vndf(40.) * r * r));
    }
    void setKernel(const KernelType kn,
                   const PS::F64 ndim) {
        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "set 1D WendlandC4!\n");
        }
    };
    
#endif

#else

#if 1
///// for kernel size independent of the number of particles (FROM)
    const PS::F64 eta     = 1.6;
//    const PS::F64 eta     = 1.6 * pow(2., 1./3.); // for r016k
//    const PS::F64 eta     = 1.6 * pow(2., 2./3.); // for r032k
//    const PS::F64 eta     = 1.6 * 2.0; // for r064k
///// for kernel size independent of the number of particles (TO)
    const PS::F64 ksrh    = 1.936492;
    const PS::F64 ksrhinv = 1. / ksrh;
    const PS::F64 dim     = 3.;
    const PS::F64 ceff0   = +3.342253804929802286e+00;
    const PS::F64 ceff1   = +1.336901521971920914e+01;

#ifdef USE_INTRINSICS

    inline vndf kernel0th(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin * rmin;
        return vndf(ceff0) * rmin2 * rmin2 * (vndf(1.) + vndf(4.) * r);
    }
    inline vndf kernel1st(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin  * rmin;
        vndf rmin3 = rmin  * rmin2;
        vndf rmin4 = rmin2 * rmin2;
        return vndf(ceff1) * (rmin4 - rmin3 * (vndf(1.) + vndf(4.) * r));
    }

#else

    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = ((1. - r > 0.) ? 1. - r : 0.);
        PS::F64 rmin2 = rmin * rmin;
        return ceff0 * rmin2 * rmin2 * (1. + 4. * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = ((1. - r > 0.) ? 1. - r : 0.);
        PS::F64 rmin2 = rmin  * rmin;
        PS::F64 rmin3 = rmin  * rmin2;
        PS::F64 rmin4 = rmin2 * rmin2;
        return ceff1 * (rmin4 - rmin3 * (1. + 4. * r));
    }

#endif

    void setKernel(const KernelType kn,
                   const PS::F64 ndim) {
        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "set 3D WendlandC2!\n");
        }
    };
#else
    const PS::F64 eta     = 1.2;
    const PS::F64 ksrh    = 2.0;
    const PS::F64 ksrhinv = 1. / ksrh;
    const PS::F64 dim     = 3.;
    const PS::F64 ceff0   = +2.5464790894703255;
    const PS::F64 ceff0x  = ceff0 * 2.;
    const PS::F64 ceff1   = -15.278874536821953;


    inline vndf kernel0th(const vndf q) {
        vndf t1 = vndf(1.)  - q;
        vndf t2 = vndf(0.5) - q;
        t1 = vndf::max(t1, 0.);
        t2 = vndf::max(t2, 0.);
        vndf t13 = t1 * t1 * t1;
        vndf t23 = t2 * t2 * t2;
        vndf w0  = t13 - vndf(4.) * t23;
        w0 *= vndf(ceff0x);        
        return w0;
    }
    inline vndf kernel1st(const vndf q) {
        vndf t1 = vndf(1.)  - q;
        vndf t2 = vndf(0.5) - q;
        vndf t3 = vndf(0.33333333333333333) - q;
        t1 = vndf::max(t1, 0.);
        t2 = vndf::max(t2, 0.);
        t3 = vndf::max(t3, 0.);
        vndf t12 = t1 * t1;
        vndf t22 = t2 * t2;
        vndf t32 = t3 * t3;
        vndf w0  = vndf::nmadd(t12, vndf(4.), t22);
        w0  = vndf::madd(w0, vndf(3.), t32);
        w0 *= vndf(ceff1);
        return w0;
    }
    void setKernel(const KernelType kn,
                   const PS::F64 ndim) {
        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "set 3D CubicSpline!\n");
        }
    };
#endif

#endif

#else
#error We have only two options: USE_IDEAL and USE_HELMHOLTZ.
#endif
    
}

namespace SK = SmoothingKernel;

/*
namespace WendlandC4 {

    const  PS::F64 eta   = 1.6;

#ifdef USE_AT1D

    const  PS::F64 dim   = +1.0;
    const  PS::F64 ksrh  = +1.936492;
    const  PS::F64 ksrhinv = 1. / ksrh;
    const  PS::F64 ceff0 = +1.5;
    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.) ? (1. - r) : 0.;
        PS::F64 rmin2 = rmin * rmin;
        return ceff0 * rmin2 * rmin2 * rmin * (1. + 5. * r + 8. * r * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.) ? (1. - r) : 0.;
        PS::F64 rmin2 = rmin  * rmin;
        PS::F64 rmin4 = rmin2 * rmin2;
        return ceff0 * rmin4 * (- (5. + 25. * r + 40. * r * r)
                                + rmin * (5. + 16. * r));
    }
    inline vndf kernel0th(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin * rmin;
        return vndf(ceff0) * rmin2 * rmin2 * rmin
            * (vndf(1.) + vndf(5.) * r + vndf(8.) * r * r);
    }
    inline vndf kernel1st(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin  * rmin;
        vndf rmin4 = rmin2 * rmin2;
        return vndf(ceff0) * rmin4 * (rmin * (vndf(5.) + vndf(16.) * r)
                                      - (vndf(5.) + vndf(25.) * r + vndf(40.) * r * r));
    }

#else

#ifdef USE_AT2D
    const  PS::F64 dim  = +2.0;
    const  PS::F64 ksrh = 2.171239;
    const  PS::F64 ksrhinv = 1. / ksrh;
    const  PS::F64 ceff0 = +2.8647889756541162e+00;
#else
    const  PS::F64 dim  = +3.0;
    const  PS::F64 ksrh = 2.207940;
    const  PS::F64 ksrhinv = 1. / ksrh;
    const  PS::F64 ceff0 = +4.9238560519055123e+00;
#endif
    const  PS::F64 ceff1 = +1.1666666666666667e+01;
    const  PS::F64 ceff2 = +2.3333333333333333e+01;
    
    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.) ? (1. - r) : 0.;
        PS::F64 rmin2 = rmin * rmin;
        return ceff0 * rmin2 * rmin2 * rmin2 * (1. + 6. * r + ceff1 * r * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.) ? (1. - r) : 0.;
        PS::F64 rmin2 = rmin  * rmin;
        PS::F64 rmin3 = rmin  * rmin2;
        return ceff0 * rmin2 * rmin3 * (- (6. + 36. * r + 70. * r * r)
                                        + rmin * (6. + ceff2 * r));
    }
    inline vndf kernel0th(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin * rmin;
        return vndf(ceff0) * rmin2 * rmin2 * rmin2
            * (vndf(1.) + vndf(6.) * r + vndf(ceff1) * r * r);
    }
    inline vndf kernel1st(const vndf r) {
        vndf rmin  = vndf::max(vndf(1.) - r, 0.);
        vndf rmin2 = rmin  * rmin;
        vndf rmin3 = rmin  * rmin2;
        return vndf(ceff0) * rmin2 * rmin3 * (rmin * (vndf(6.) + vndf(ceff2) * r)
                                              - (vndf(6.) + vndf(36.) * r + vndf(70.) * r * r));
    }
    
#endif
    
}
*/
