# 800/200 = 4
BEGIN{
    n = 4;
}
{
    if(NR == 1) {
        print $4, $5, $6;
    } else {
        dz = ($4 - prevz) / n;
        dd = ($5 - prevd) / n;
        dv = ($6 - prevv) / n;
        for(i = 0; i < n; i++) {
            posz = prevz + dz * (i + 1);
            dens = prevd + dd * (i + 1);
            velz = prevv + dv * (i + 1);
            printf("%+e %+e %+e\n", posz, dens, velz);
        }
    }
    prevz = $4;
    prevd = $5;
    prevv = $6;
}
END{
    for(i = 0; i < n-1; i++) {
        posz = prevz + dz * (i + 1);
        dens = prevd + dd * (i + 1);
        velz = prevv + dv * (i + 1);
        printf("%+e %+e %+e\n", posz, dens, velz);
    }
}

