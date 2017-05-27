BEGIN{
    pi  = 3.141;
    idc = 0;
}
{
    id    = $1;
    istar = $2;
    mass  = $3;
    px    = $4;
    py    = $5;
    pz    = $6;
    vx    = $7;
    vy    = $8;
    vz    = $9;
    ue    = $13;
    alph  = $14;
    alphu = $15;
    ksr   = $17;
    for(k = 0; k < 13; k++) {
        cmps[k] = $(32+k);
    }

    phi = pi * (2. * rand() - 1.);
    cth = 2. * rand() - 1.;
    sth = sqrt(1. - cth**2);
    dr  = (1./2.)**(1/3.) * (ksr * 0.1);
    dx  = dr * sth * cos(phi);
    dy  = dr * sth * sin(phi);
    dz  = dr * cth;
    dvx = 0. * sth * cos(phi);
    dvy = 0. * sth * sin(phi);
    dvz = 0. * cth;

    id0    = idc;
    istar0 = istar;
    mass0  = 0.5 * mass;
    px0    = px + dx;
    py0    = py + dy;
    pz0    = pz + dz;
    vx0    = vx + dvx;
    vy0    = vy + dvy;
    vz0    = vz + dvz;
    ue0    = ue;
    alph0  = alph;
    alphu0 = alphu;
    ksr0   = ksr;
    for(k = 0; k < 13; k++) {
        cmps0[k] = cmps[k]
    }
    idc++;

    id1    = idc;
    istar1 = istar;
    mass1  = 0.5 * mass;
    px1    = px - dx;
    py1    = py - dy;
    pz1    = pz - dz;
    vx1    = vx - dvx;
    vy1    = vy - dvy;
    vz1    = vz - dvz;
    ue1    = ue;
    alph1  = alph;
    alphu1 = alphu;
    ksr1   = ksr;
    for(k = 0; k < 13; k++) {
        cmps1[k] = cmps[k]
    }
    idc++;

    printf("%10d %2d %+e", id0, istar0, mass0);
    printf(" %+e %+e %+e", px0, py0, pz0);
    printf(" %+e %+e %+e", vx0, vy0, vz0);
    printf(" %+e %+e %+e", ue0, alph0, alphu0);
    printf(" %+e", ksr0);
    for(k = 0; k < 13; k++) {
        printf(" %+.3e", cmps0[k]);
    }
    printf("\n");

    printf("%10d %2d %+e", id1, istar1, mass1);
    printf(" %+e %+e %+e", px1, py1, pz1);
    printf(" %+e %+e %+e", vx1, vy1, vz1);
    printf(" %+e %+e %+e", ue1, alph1, alphu1);
    printf(" %+e", ksr1);
    for(k = 0; k < 13; k++) {
        printf(" %+.3e", cmps1[k]);
    }
    printf("\n");
}
