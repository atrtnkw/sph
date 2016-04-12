BEGIN{
    sum = 0.;
}
{
    sum += $3 * (0.5 * ($7**2 + $8**2 + $9**2) + $45);
}
END{
    printf("%+e\n", sum);
}


