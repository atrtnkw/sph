BEGIN{
    rad0 = 0.;
    menc = 0.;
    msun = 1.989e33;
}
{
    rad1 = $1;
    dens = $2;
    temp = $3;
    eint = $4;
    pres = $5;
    velx = $6;

    menc += 4. * 3.141 / 3. * (rad1**3 - rad0**3) * dens

    printf("%+e %+e %+e %+e %+e %+e %+e", rad1, dens, pres, temp, eint, menc/msun, velx);
    for(k = 0; k < 13; k++) {
        printf(" %+e", $(8+k));
    }
    printf("\n");

    rad0 = rad1;
}
