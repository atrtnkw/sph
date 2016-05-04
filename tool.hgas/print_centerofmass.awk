BEGIN {
    px = 0.;
    py = 0.;
    pz = 0.;
    vx = 0.;
    vy = 0.;
    vz = 0.;
    np = 0;
}
{
    px += $4;
    py += $5;
    pz += $6;
    vx += $7;
    vy += $8;
    vz += $9;
    np += 1;
}
END {
    ninv = 1. / np;
    px *= ninv;
    py *= ninv;
    pz *= ninv;
    vx *= ninv;
    vy *= ninv;
    vz *= ninv;
    printf("x: %+e %+e %+e\n", px, py, pz);
    printf("v: %+e %+e %+e\n", vx, vy, vz);
}
