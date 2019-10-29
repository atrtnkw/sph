{
    if(NR % pitch == 0) {
        printf("%+e %+e %+e\n", $1, $2, -$3);
    }
}
