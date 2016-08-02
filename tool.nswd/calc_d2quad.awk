BEGIN{
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            ii = 3 * i + j;
            q[ii] = 0.;
        }
    }
}
{
    m    = $3;
    x[0] = $4  - cx;
    x[1] = $5  - cy;
    x[2] = $6  - cz;
    v[0] = $7  - vx;
    v[1] = $8  - vy;
    v[2] = $9  - vz;
    a[0] = $10 - ax;
    a[1] = $11 - ay;
    a[2] = $12 - az;
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            if(j < i) {
                continue;
            }
            qq = 2. * v[i] * v[j] + x[i] * a[j] +  a[i] * x[j];
            ii = 3 * i + j;
            q[ii] += m * qq;
        }
    }
}
END{
    c = 2.998e10;
    g = 6.673e-8;
    d = 3.086e18 * 1e4;
    ceff = g/(c**4*d);
    printf("%16.10f", time);
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            if(j < i) {
                continue;
            }
            ii = 3 * i + j;
#            printf(" %+.16e", q[ii]);            
            printf(" %+.16e", q[ii] * g / c**4);
        }
    }
#    printf(" %+e %+e", ceff*(q[0]-q[4]), ceff*(2.*q[1]));
    printf("\n");
}
