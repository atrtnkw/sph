{
    printf("%8d %2d %+e", $1, $2, $3);
    printf(" %+e %+e %+e", $4, $5, $6);
    printf(" %+e %+e %+e", $7, $8, $9);
    printf(" %+e %+e %+e", $10, $11, 0.);
    printf(" %+e", $13);
    for(i = 0; i < 13; i++) {
        printf(" %+.3e", $(i+14));
    }
    printf("\n");
}
