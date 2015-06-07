#pragma once

namespace WendlandC2 {

    const  PS::F64 eta   = 1.6d;

#ifdef USE_AT1D

    const  PS::F64 dim   = +1.0d;
    const  PS::F64 ksrh  = 1.620185d;
    const  PS::F64 ceff0 = +1.25d;
    const  PS::F64 ceff1 = +3.75d;
    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.d) ? (1.d - r) : 0.;
        return ceff0 * rmin * rmin * rmin * (1.d + 3.d * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.) ? (1. - r) : 0.;
        PS::F64 rmin2 = rmin  * rmin;
        PS::F64 rmin3 = rmin  * rmin2;
        return ceff1 * (rmin3 - rmin2 * (1.d + 3.d * r));
    }

#else

#ifdef USE_AT2D
    const  PS::F64 dim  = +2.0;
    const  PS::F64 ksrh = 1.897367d;
    const  PS::F64 ceff0 = +2.228169203286535005e+00d;
    const  PS::F64 ceff1 = +8.912676813146140020e+00d;
#else
    const  PS::F64 dim  = +3.0;
    const  PS::F64 ksrh = 1.936492d;
    const  PS::F64 ceff0 = +3.342253804929802286e+00d;
    const  PS::F64 ceff1 = +1.336901521971920914e+01d;
#endif
    
    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.d) ? (1.d - r) : 0.;
        PS::F64 rmin2 = rmin * rmin;
        return ceff0 * rmin2 * rmin2 * (1.d + 4.d * r);
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 rmin  = (r < 1.) ? (1. - r) : 0.;
        PS::F64 rmin2 = rmin  * rmin;
        PS::F64 rmin3 = rmin  * rmin2;
        PS::F64 rmin4 = rmin2 * rmin2;        
        return ceff1 * (rmin4 - rmin3 * (1.d + 4.d * r));
    }
    
#endif
    
}

namespace WendlandC4 {

    const  PS::F64 eta   = 1.6d;

#ifdef USE_AT1D

    const  PS::F64 dim   = +1.0d;
    const  PS::F64 ksrh  = +1.936492d;
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

#else

#ifdef USE_AT2D
    const  PS::F64 dim  = +2.0;
    const  PS::F64 ksrh = 2.171239d;
    const  PS::F64 ceff0 = +2.8647889756541162e+00d;
#else
    const  PS::F64 dim  = +3.0;
    const  PS::F64 ksrh = 2.207940d;
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
    
#endif
    
}

namespace CubicSpline {
    
    const PS::F64 eta   = 1.2d;
//    const PS::F64 eta   = 4.8d;

#ifdef USE_AT1D
    const PS::F64 dim   = +1.0;
    const PS::F64 ceff0 = +1.3333333333333333d;
    const PS::F64 ceff1 = -8.d;
#elif defined USE_AT2D
    const PS::F64 dim   = +2.0;
    const PS::F64 ceff0 = +1.8189136353359468d;
    const PS::F64 ceff1 = -10.913481812015681d;
#else
    const PS::F64 dim   = +3.0;
    const PS::F64 ceff0 = +2.5464790894703255d;
    const PS::F64 ceff1 = -15.278874536821953d;
#endif

    const  PS::F64 ksrh = 2.0d;
    
    inline PS::F64 kernel0th(const PS::F64 r) {
        PS::F64 w0 = 0.d;
        if(r < 0.5d) {
            w0 = 1.d - 6.d * r * r + 6.d * r * r * r;
        } else if(r < 1.d) {
            w0 = 2.d * (1.d - r) * (1.d - r) * (1.d - r);
        }
        return ceff0 * w0;
    }
    inline PS::F64 kernel1st(const PS::F64 r) {
        PS::F64 w1 = 0.d;
        if(r < 0.33333333333333333d) {
            w1 = 0.33333333333333333d;
        } else if(r < 0.5d) {
            w1 = 2.d * r - 3.d * r * r;
        } else if(r < 1.0d) {
            w1 = (1.d - r) * (1.d - r);
        }
        return ceff1 * w1;
    }
    
}

#ifdef WENDLANDC2
namespace KernelSph = WendlandC2;
#elif defined WENDLANDC4
namespace KernelSph = WendlandC4;
#elif defined CUBICSPLINE
namespace KernelSph = CubicSpline;
#else
#error
#endif
