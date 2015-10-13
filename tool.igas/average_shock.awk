#argv: xmin, width
BEGIN{
    n = -1;
}
{
    i = int(($3 - xmin) / width);
    nbin[i]   += 1;
    vxbin[i]  += $6;
    rhobin[i] += $16;
    pbin[i]   += $18;
    csbin[i]  += $17;
    if(i > n) {
        n = i;
    }
}
END{
    for(i = 0; i < n; i++) {
        vx  = vxbin[i]  / nbin[i];
        rho = rhobin[i] / nbin[i];
        p   = pbin[i]    / nbin[i];
        cs  = csbin[i]  / nbin[i]  
            printf("%+e %+e %+e %+e %+e\n", width * (i + 0.5) + xmin, vx, rho, p, cs);
    }
}
