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
    printf(" %+e", mele[0]/msun);
    printf(" %+e", mele[1]/msun);
    printf(" %+e", mele[2]/msun);
    printf(" %+e", mele[3]/msun);
    printf(" %+e", mele[4]/msun);
    printf(" %+e", (mele[6]+mele[7]+mele[8]+mele[9]+mele[10])/msun);
    printf(" %+e", (mele[11]+mele[12])/msun);
    printf("\n");
}
