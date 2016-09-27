# argv: n1msun
BEGIN{
    pi = 3.14159265359;
    ms = 1.989e33;
    
    cdz  = 1e-2;
    rin  = 1e6;
    rout = 1e11;
    ndig = 10
    dr   = 10**(1./ndig);
    n    = int(10 * log(rout / rin)/log(10.));
    for(i = 0; i < n; i++) {
	sigma[i] = 0.0;
	rhocn[i] = 0.0;
	tempn[i] = 0.0;
	ncn[i]   = 0;
    }

}
{
    r = sqrt($4**2 + $5**2);
    if(r < rin || rout < r) {
	next;
    }    
    i = int(log(r/rin) / log(dr));
    sigma[i]++;
    ri = rin * dr**i;
    if(-cdz*ri < $6 && $6 < cdz*ri) {
	rhocn[i] += $16;
	tempn[i] += $21;
	ncn[i]++;
    }
}
END{
    for(i = 0; i < n; i++) {
	ri = rin * dr**i;
	rn = rin * dr**(i+1);
	si = sigma[i] * ms / n1msun / (4. * pi * (rn**2 - ri**2));
	ni = ((ncn[i] != 0) ? 1 / ncn[i] : 0.);
	dc = rhocn[i] * ni;
	tc = tempn[i] * ni;
	ht = ((dc != 0.) ? 0.5 * si / (dc * ri) : 0.);
	yi = ((ni != 0 && tc > 1e8) ? calcYcDotInv(dc, tc) : 0.);
	printf("%+e %+e %+e %+e %+e %+e %9d\n", ri, si, dc, tc, ht, yi, sigma[i]);
    }
}

function calcYcDotInv(rho, temp)
{
    qc = 4.48e18;
    yc = 0.033;
    lq = 84.165;

    t9  = temp * 1e-9;
    t9a = t9 / (1. + 0.067 * t9);
    tAA = 8.54e26 * t9a**(5./6.) / (t9**1.5);

    ycdot = 1. / (rho * qc * tAA * yc * yc * exp(-lq / t9a**(1./3.)));

    return ycdot;
}


