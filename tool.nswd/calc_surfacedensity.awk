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
	printf("%+e %+e %+e %9d\n", ri, si, dc, sigma[i]);
    }
}

