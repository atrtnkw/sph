#pragma once

#include <immintrin.h>

struct v8sf {
    typedef float _v8sf __attribute__((vector_size(32))) __attribute__((aligned(32)));
    _v8sf val;

    static int getVectorLength() {
        return 8;
    }

    v8sf(const v8sf & rhs) : val(rhs.val) {}
    v8sf operator = (const v8sf rhs) {
        val = rhs.val;
        return (*this);
    }

    v8sf() : val(_mm256_setzero_ps()) {}
    v8sf(const float x) : val(_mm256_set_ps(x, x, x, x, x, x, x, x)) {}
    v8sf(const float x0, const float y0, const float z0, const float w0,
         const float x1, const float y1, const float z1, const float w1)
        : val(_mm256_set_ps(w1, z1, y1, x1, w0, z0, y0, x0)) {}

    v8sf(const _v8sf _val) : val(_val) {}
    operator _v8sf() {return val;}

    v8sf operator + (const v8sf rhs) const {
        return v8sf(_mm256_add_ps(val, rhs.val));
    }
    v8sf operator - (const v8sf rhs) const {
        return v8sf(_mm256_sub_ps(val, rhs.val));
    }
    v8sf operator * (const v8sf rhs) const {
        return v8sf(_mm256_mul_ps(val, rhs.val));
    }
    v8sf operator / (const v8sf rhs) const {
        return v8sf(_mm256_div_ps(val, rhs.val));
    }
    static v8sf madd(const v8sf c, const v8sf a, const v8sf b) {
        return v8sf(_mm256_fmadd_ps(a.val, b.val, c.val));
    }
    static v8sf nmadd(const v8sf c, const v8sf a, const v8sf b) {
        return v8sf(_mm256_fnmadd_ps(a.val, b.val, c.val));
    }

    v8sf operator += (const v8sf rhs) {
        val = _mm256_add_ps(val, rhs.val);
        return (*this);
    }
    v8sf operator -= (const v8sf rhs) {
        val = _mm256_sub_ps(val, rhs.val);
        return (*this);
    }
    v8sf operator *= (const v8sf rhs) {
        val = _mm256_mul_ps(val, rhs.val);
        return (*this);
    }
    v8sf operator /= (const v8sf rhs) {
        val = _mm256_div_ps(val, rhs.val);
        return (*this);
    }

    static v8sf sqrt(const v8sf rhs) {
        return v8sf(_mm256_sqrt_ps(rhs.val));
    }

    static v8sf rsqrt(const v8sf rhs) {
        return v8sf(_mm256_rsqrt_ps(rhs.val));
    }

    void store(float *p) const {
        _mm256_store_ps(p, val);
    }

    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n%+e %+e %+e %+e\n\n") const {
        int vl = getVectorLength();
        float a[vl];
        store(a);
        fprintf(fp, fmt, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);
    }    
};

struct v4df {
    typedef double _v4df __attribute__((vector_size(32))) __attribute__((aligned(32)));
    _v4df val;

    static int getVectorLength() {
        return 4;
    }

    v4df(const v4df & rhs) : val(rhs.val) {}
    v4df operator = (const v4df rhs) {
        val = rhs.val;
        return (*this);
    }

    v4df() : val(_mm256_setzero_pd()) {}
    v4df(const double x) : val(_mm256_set_pd(x, x, x, x)) {}
    v4df(const double x, const double y, const double z, const double w)
        : val(_mm256_set_pd(w, z, y, x)) {}

    v4df(const _v4df _val) : val(_val) {}
    operator _v4df() {return val;}

    v4df operator + (const v4df rhs) const {
        return v4df(_mm256_add_pd(val, rhs.val));
    }
    v4df operator - (const v4df rhs) const {
        return v4df(_mm256_sub_pd(val, rhs.val));
    }
    v4df operator * (const v4df rhs) const {
        return v4df(_mm256_mul_pd(val, rhs.val));
    }
    v4df operator / (const v4df rhs) const {
        return v4df(_mm256_div_pd(val, rhs.val));
    }

    static v4df madd(const v4df c, const v4df a, const v4df b) {
        return v4df(_mm256_fmadd_pd(a.val, b.val, c.val));
    }
    static v4df nmadd(const v4df c, const v4df a, const v4df b) {
        return v4df(_mm256_fnmadd_pd(a.val, b.val, c.val));
    }
    static v4df max(const v4df a, const v4df b) {
        return v4df(_mm256_max_pd(a.val, b.val));
    }
    static v4df min(const v4df a, const v4df b) {
        return v4df(_mm256_min_pd(a.val, b.val));
    }

    v4df operator += (const v4df rhs) {
        val = _mm256_add_pd(val, rhs.val);
        return (*this);
    }
    v4df operator -= (const v4df rhs) {
        val = _mm256_sub_pd(val, rhs.val);
        return (*this);
    }
    v4df operator *= (const v4df rhs) {
        val = _mm256_mul_pd(val, rhs.val);
        return (*this);
    }
    v4df operator /= (const v4df rhs) {
        val = _mm256_div_pd(val, rhs.val);
        return (*this);
    }

    v4df operator & (const v4df rhs) {
        return v4df(_mm256_and_pd(val, rhs.val));
    }
    v4df operator != (const v4df rhs) {
        return v4df(_mm256_cmp_pd(val, rhs.val, _CMP_NEQ_UQ));
    }
    v4df operator < (const v4df rhs) {
        return v4df(_mm256_cmp_pd(val, rhs.val, _CMP_LT_OS));
    }

    static v4df sqrt(const v4df rhs) {
        return v4df(_mm256_sqrt_pd(rhs.val));
    }

    void store(double *p) const {
        _mm256_store_pd(p, val);
    }
    void load(double const *p) {
        val = _mm256_load_pd(p);
    }

    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n") const {
        int vl = getVectorLength();
        double a[vl];
        _mm256_store_pd(a, val);
        fprintf(fp, fmt, a[0], a[1], a[2], a[3]);
    }    
};

