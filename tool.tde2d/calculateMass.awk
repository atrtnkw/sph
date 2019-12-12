BEGIN {
    for(i = 0; i < 14; i++) {
	mass[i] = 0.;
    }
}
{
    mass[0] += $5;
    for(i = 1; i < 14; i++) {
	mass[i] += $5 * $(i+6)
    }
}
END {
    for(i = 0; i < 14; i++) {
	printf("%2d %+e\n", i, mass[i]);
    }
}
