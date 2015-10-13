BEGIN{
    n = 0;
}
{
    if(NR > 1) {
        id[n] = $1;
        is[n] = $2;
        ms[n] = $3;
        px[n] = $4;
        py[n] = $5;
        pz[n] = $6;
        vx[n] = $7;
        vy[n] = $8;
        vz[n] = $9;
        eg[n] = $10;
        al[n] = $11;
        hs[n] = $12;
        n++;
    }
}
END{
    time = 0.0;
    tend = 100.0;
    dtsp = 1.0;
    alphamax = 2.0;
    alphamin = 0.1;
    dummy    = 3e6;
    tceff    = 0.1;
        
    printf("%+e %+e %+e\n", time, tend, dtsp);
    printf("%+e %+e %+e\n", alphamax, alphamin, tceff);
    printf("%+e\n", dummy);
    printf("%8d\n", n);
    for(i = 0; i < n; i++) {
        printf("%8d %2d %+.16e", id[i], is[i], ms[i]);
        printf(" %+.16e %+.16e %+.16e", px[i], py[i], pz[i]);
        printf(" %+.16e %+.16e %+.16e", vx[i], vy[i], vz[i]);
        printf(" %+.16e %+.16e %+.16e %+.16e", eg[i], al[i], al[i], hs[i]);
        printf("\n");
    }
}

