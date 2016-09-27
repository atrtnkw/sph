BEGIN{
    dndex = 8;
    nfbin = 14 * dndex + 1;
    pi = 3.14159265359;
    for(i = 0; i < nfbin; i++) {
	x[i] = 0.0;
	y[i] = 0.0;
	f[i] = 2**(i/dndex-7);
#	printf("%16.10f\n", f[i]);
    }
    n = 0;
}
{
    if($1 < tbgn || tend < $1) {
	next;
    }

    xt    = $2 - $5;
    for(i = 0; i < nfbin; i++) {
	x[i] += xt * cos(- 2.0 * pi * f[i] * $1);
	y[i] += xt * sin(- 2.0 * pi * f[i] * $1);
    }
    n++;
}
END{
    for(i = 0; i < nfbin; i++) {
	x[i] /= n;
	y[i] /= n;
	p = sqrt(x[i]**2 + y[i]**2);
	printf("%+e %+e\n", f[i], p);
    }
}
