# 800/200 = 4
BEGIN{
    n = 4;
}
{
    if(NR == 1) {
        print $1, $2, $3;
    } else {
        dz = ($1 - prevz) / n;
        dd = ($2 - prevd) / n;
        dv = ($3 - prevv) / n;
        for(i = 0; i < n; i++) {
            posz = prevz + dz * (i + 1);
            dens = prevd + dd * (i + 1);
            velz = prevv + dv * (i + 1);
            printf("%+e %+e %+e\n", posz, dens, velz);
        }
    }
    prevz = $1;
    prevd = $2;
    prevv = $3;
}
END{
    for(i = 0; i < n-1; i++) {
        posz = prevz + dz * (i + 1);
        dens = prevd + dd * (i + 1);
        velz = prevv + dv * (i + 1);
        printf("%+e %+e %+e\n", posz, dens, velz);
    }
}

