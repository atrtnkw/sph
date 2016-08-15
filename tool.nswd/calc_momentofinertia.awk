BEGIN{
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            ii = 3 * i + j;
            q[ii] = 0.;
        }
    }
}
{
    m    = $3;
    x[0] = $4  - cx;
    x[1] = $5  - cy;
    x[2] = $6  - cz;
    x2   = x[0]**2 + x[1]**2 + x[2]**2;
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            if(j < i) {
                continue;
            }
            qq = - x[i] * x[j] + ((j == i) ? x2 : 0.);
            ii = 3 * i + j;
            q[ii] += m * qq;
        }
    }
}
END{
    printf("%16.10f", time);
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            if(j < i) {
                continue;
            }
            ii = 3 * i + j;
            printf(" %+.16e", q[ii]);            
        }
    }
    printf("\n");
}
