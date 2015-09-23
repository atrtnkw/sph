# argv: id (0/1), time
BEGIN{
    xc[0] = xc[1] = xc[2] = 0.0;
    vc[0] = vc[1] = vc[2] = 0.0;
    n  = 0;
}
{
    if(NR == 1) {
        if(id == 0) {
            printf("Output primary angz\n");
        } else {
            printf("Output secondary angz\n");
        }
    }

    if($2 != id)
        next;

    px[n]  = $4;
    py[n]  = $5;
    pz[n]  = $6;
    vx[n]  = $7;
    vy[n]  = $8;
    vz[n]  = $9;

    xc[0] += px[n];
    xc[1] += py[n];
    xc[2] += pz[n];
    vc[0] += vx[n];
    vc[1] += vy[n];
    vc[2] += vz[n];

    n++;
}
END{
    ninv   = 1. / n;
    xc[0] *= ninv;
    xc[1] *= ninv;
    xc[2] *= ninv;
    vc[0] *= ninv;
    vc[1] *= ninv;
    vc[2] *= ninv;

    ang[0] = ang[1] = ang[2] = 0.0;
    for(i = 0; i < n; i++) {
        dpx  = px[i] - xc[0];
        dpy  = py[i] - xc[1];
        dpz  = pz[i] - xc[2];
        dvx  = vx[i] - vc[0];
        dvy  = vy[i] - vc[1];
        dvz  = vz[i] - vc[2];
        ang[0] += dpy * dvz - dpz * dvy;
        ang[1] += dpz * dvx - dpx * dvz;
        ang[2] += dpx * dvy - dpy * dvx;
    }
    printf("%4d %+e %+e %+e\n", time, ang[0], ang[1], ang[2]);
}


