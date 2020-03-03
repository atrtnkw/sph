BEGIN {
    name[0]  = "Al";
    name[1]  = "He";
    name[2]  = "C";
    name[3]  = "O";
    name[4]  = "Ne";
    name[5]  = "Mg";
    name[6]  = "Si";
    name[7]  = "S";
    name[8]  = "Ar";
    name[9]  = "Ca";
    name[10] = "Ti";
    name[11] = "Cr";
    name[12] = "Fe";
    name[13] = "Ni";
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
        printf("%2s %+e %+e\n", name[i], mass[i], mass[i]/mass[0]);
    }
}
