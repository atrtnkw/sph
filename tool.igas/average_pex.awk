#argv: width
BEGIN{
    n = -1;
}
{
    r  = sqrt($3**2+$4**2+$5**2);
    vr = ($3 * $6 + $4 * $7 + $5 * $8) / r;
    i = int(r / width);
    nbin[i]   += 1;
    vxbin[i]  += vr;
    rhobin[i] += $15;
    pbin[i]   += $17;
    csbin[i]  += $16;
    if(i > n) {
        n = i;
    }
}
END{
    for(i = 0; i < n; i++) {
        if(nbin[i] == 0)
            continue
        vx  = vxbin[i]  / nbin[i];
        rho = rhobin[i] / nbin[i];
        p   = pbin[i]    / nbin[i];
        cs  = csbin[i]  / nbin[i]  
            printf("%+e %+e %+e %+e %+e\n", width * (i + 0.5) + xmin, vx, rho, p, cs);
    }
}
