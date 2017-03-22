#argv: dtsp
BEGIN{
    msun = 1.989e33;
    nele = 13;
    for(i = 0; i < nele; i++) {
        mele[i] = 0.;
    }
}
{
    for(i = 0; i < nele; i++) {
        mele[i] += $3 * $(32+i);
    }    
}
END{
    printf("%16.8f", time/dtsp);
    for(i = 0; i < nele; i++) {
        printf(" %+e", mele[i]/msun);
    }
    printf("\n");
}
