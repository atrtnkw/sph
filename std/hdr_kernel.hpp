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

    v4df (*kernel0th)(const v4df q);
    v4df (*kernel1st)(const v4df q);

    inline v4df calcCubicSpline0th(const v4df q) {
        v4df t1 = v4df(1.)  - q;
        v4df t2 = v4df(0.5) - q;
        t1 = v4df::max(t1, 0.);
        t2 = v4df::max(t2, 0.);
        v4df t13 = t1 * t1 * t1;
        v4df t23 = t2 * t2 * t2;
        v4df w0  = t13 - v4df(4.) * t23;
        w0 *= v4df(ceff0x);
        
        return w0;
    }
    inline v4df calcCubicSpline1st(const v4df q) {
        v4df t1 = v4df(1.)  - q;
        v4df t2 = v4df(0.5) - q;
        v4df t3 = v4df(0.33333333333333333) - q;
        t1 = v4df::max(t1, 0.);
        t2 = v4df::max(t2, 0.);
        t3 = v4df::max(t3, 0.);
        v4df t12 = t1 * t1;
        v4df t22 = t2 * t2;
        v4df t32 = t3 * t3;
        v4df w0  = v4df::nmadd(t12, v4df(4.), t22);
        w0  = v4df::madd(w0, v4df(3.), t32);
        w0 *= v4df(ceff1);
        return w0;
    }

    inline v4df calcWendlandC20th1D(const v4df r) {
        v4df rmin = v4df::max(v4df(1.d) - r, 0.d);
        return v4df(ceff0) * rmin * rmin * rmin * (v4df(1.d) + v4df(3.d) * r);
    }
    inline v4df calcWendlandC21st1D(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin * rmin;
        v4df rmin3 = rmin * rmin2;
        return v4df(ceff1) * (rmin3 - rmin2 * (v4df(1.d) + v4df(3.d) * r));
    }

    inline v4df calcWendlandC20thMD(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin * rmin;
        return v4df(ceff0) * rmin2 * rmin2 * (v4df(1.d) + v4df(4.d) * r);
    }
    inline v4df calcWendlandC21stMD(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin  * rmin;
        v4df rmin3 = rmin  * rmin2;
        v4df rmin4 = rmin2 * rmin2;
        return v4df(ceff1) * (rmin4 - rmin3 * (v4df(1.d) + v4df(4.d) * r));
    }

    void setKernel(const KernelType kn,
                   const PS::F64 ndim) {
        if(kn == CubicSpline) {
            eta       = 1.2;
            ksrh      = 2.0d;
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
                ceff0 = +1.8189136353359468d;
                ceff1 = -10.913481812015681d;
                if(PS::Comm::getRank() == 0) {
                    fprintf(stderr, "set 2D cubic spline!\n");
                }
            } else {
                dim   = +3.0;
                ceff0  = +2.5464790894703255d;
                ceff1  = -15.278874536821953d;
                if(PS::Comm::getRank() == 0) {
                    fprintf(stderr, "set 3D cubic spline!\n");
                }
            }
            ceff0x = ceff0 * 2.d;
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

    const PS::F64 eta     = 1.6;
    const PS::F64 ksrh    = 1.620185;
    const PS::F64 ksrhinv = 1. / ksrh;
    const PS::F64 dim     = 1.;
    const PS::F64 ceff0   = +1.25;
    const PS::F64 ceff1   = +3.75;

    inline v4df kernel0th(const v4df r) {
        v4df rmin = v4df::max(v4df(1.d) - r, 0.d);
        return v4df(ceff0) * rmin * rmin * rmin * (v4df(1.d) + v4df(3.d) * r);
    }
    inline v4df kernel1st(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin * rmin;
        v4df rmin3 = rmin * rmin2;
        return v4df(ceff1) * (rmin3 - rmin2 * (v4df(1.d) + v4df(3.d) * r));
    }
    void setKernel(const KernelType kn,
                   const PS::F64 ndim) {
        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "set 1D WendlandC2!\n");
        }
    };

#else

    const PS::F64 eta     = 1.6;
    const PS::F64 ksrh    = 1.936492;
    const PS::F64 ksrhinv = 1. / ksrh;
    const PS::F64 dim     = 3.;
    const PS::F64 ceff0   = +3.342253804929802286e+00;
    const PS::F64 ceff1   = +1.336901521971920914e+01;

    inline v4df kernel0th(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin * rmin;
        return v4df(ceff0) * rmin2 * rmin2 * (v4df(1.d) + v4df(4.d) * r);
    }
    inline v4df kernel1st(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin  * rmin;
        v4df rmin3 = rmin  * rmin2;
        v4df rmin4 = rmin2 * rmin2;
        return v4df(ceff1) * (rmin4 - rmin3 * (v4df(1.d) + v4df(4.d) * r));
    }
    void setKernel(const KernelType kn,
                   const PS::F64 ndim) {
        if(PS::Comm::getRank() == 0) {
            fprintf(stderr, "set 3D WendlandC2!\n");
        }
    };

#endif

#else
#error We have only two options: USE_IDEAL and USE_HELMHOLTZ.
#endif
    
}

namespace SK = SmoothingKernel;

/*
namespace WendlandC4 {

    const  PS::F64 eta   = 1.6d;

#ifdef USE_AT1D

    const  PS::F64 dim   = +1.0d;
    const  PS::F64 ksrh  = +1.936492d;
    const  PS::F64 ksrhinv = 1. / ksrh;
    const  PS::F64 ceff0 = +1.5d;
    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.d) ? (1.d - r) : 0.d;
        PS::F64 rmin2 = rmin * rmin;
        return ceff0 * rmin2 * rmin2 * rmin * (1.d + 5.d * r + 8.d * r * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.d) ? (1.d - r) : 0.d;
        PS::F64 rmin2 = rmin  * rmin;
        PS::F64 rmin4 = rmin2 * rmin2;
        return ceff0 * rmin4 * (- (5.d + 25.d * r + 40.d * r * r)
                                + rmin * (5.d + 16.d * r));
    }
    inline v4df kernel0th(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin * rmin;
        return v4df(ceff0) * rmin2 * rmin2 * rmin
            * (v4df(1.d) + v4df(5.d) * r + v4df(8.d) * r * r);
    }
    inline v4df kernel1st(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin  * rmin;
        v4df rmin4 = rmin2 * rmin2;
        return v4df(ceff0) * rmin4 * (rmin * (v4df(5.d) + v4df(16.d) * r)
                                      - (v4df(5.d) + v4df(25.d) * r + v4df(40.d) * r * r));
    }

#else

#ifdef USE_AT2D
    const  PS::F64 dim  = +2.0;
    const  PS::F64 ksrh = 2.171239d;
    const  PS::F64 ksrhinv = 1. / ksrh;
    const  PS::F64 ceff0 = +2.8647889756541162e+00d;
#else
    const  PS::F64 dim  = +3.0;
    const  PS::F64 ksrh = 2.207940d;
    const  PS::F64 ksrhinv = 1. / ksrh;
    const  PS::F64 ceff0 = +4.9238560519055123e+00d;
#endif
    const  PS::F64 ceff1 = +1.1666666666666667e+01d;
    const  PS::F64 ceff2 = +2.3333333333333333e+01d;
    
    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.d) ? (1.d - r) : 0.;
        PS::F64 rmin2 = rmin * rmin;
        return ceff0 * rmin2 * rmin2 * rmin2 * (1.d + 6.d * r + ceff1 * r * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.) ? (1. - r) : 0.;
        PS::F64 rmin2 = rmin  * rmin;
        PS::F64 rmin3 = rmin  * rmin2;
        return ceff0 * rmin2 * rmin3 * (- (6.d + 36.d * r + 70.d * r * r)
                                        + rmin * (6.d + ceff2 * r));
    }
    inline v4df kernel0th(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin * rmin;
        return v4df(ceff0) * rmin2 * rmin2 * rmin2
            * (v4df(1.d) + v4df(6.d) * r + v4df(ceff1) * r * r);
    }
    inline v4df kernel1st(const v4df r) {
        v4df rmin  = v4df::max(v4df(1.d) - r, 0.d);
        v4df rmin2 = rmin  * rmin;
        v4df rmin3 = rmin  * rmin2;
        return v4df(ceff0) * rmin2 * rmin3 * (rmin * (v4df(6.d) + v4df(ceff2) * r)
                                              - (v4df(6.d) + v4df(36.d) * r + v4df(70.d) * r * r));
    }
    
#endif
    
}
*/
