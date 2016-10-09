BEGIN{
    a    = 2.0e-8; # density fitting
    b    = 1.2e+8;  # density fitting
    c    = 24.;    # velocity fitting
    d    = -50.;   # energy fitting
    e    = 5.1e17; # energy fitting
    n    = 1000;
    zmax = 8e7;

    mtot = -a/3. * zmax**3 + b * zmax;
    m    = mtot / n;
    for(i = 0; i < n; i++) {
        if(i == 0) {
            pz[0] = 0.;
        } else {
            zp    = pz[i-1];
            pz[i] = zp + m / (-a * zp**2 + b);
        }
        vz[i] = -c * pz[i];
        ue[i] =  d * pz[i]**2 + e;
        hs[i] = 5. * m / (-a * zp**2 + b);
    }
    for(i = n-1; i >= 0; i--) {
        if(i == 0) {
            continue;
        }
        printf(" %8d %2d %+e", n-1-i, 0, m);
        printf(" %+e %+e %+e", -pz[i], 0.0, 0.0);
        printf(" %+e %+e %+e", -vz[i], 0.0, 0.0);
        printf(" %+e %+e %+e",  ue[i], 1.0, 0.0);
        printf(" %+e", hs[i]);
        printf(" 0.0 0.6 0.35 0.05");
        for(k = 0; k < 9; k++) {
            printf(" 0.0");
        }
        printf("\n");
    }
    for(i = 0; i < n; i++) {
        printf(" %8d %2d %+e", i+n-1, 0, m);
        printf(" %+e %+e %+e", pz[i], 0.0, 0.0);
        printf(" %+e %+e %+e", vz[i], 0.0, 0.0);
        printf(" %+e %+e %+e", ue[i], 1.0, 0.0);
        printf(" %+e", hs[i]);
        printf(" 0.0 0.6 0.35 0.05");
        for(k = 0; k < 9; k++) {
            printf(" 0.0");
        }
        printf("\n");
    }
}
