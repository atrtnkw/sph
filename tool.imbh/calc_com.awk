BEGIN{
    xc = 0.;
    yc = 0.;
    zc = 0.;
    nc = 0;
}
{
    xc += $4;
    yc += $5;
    zc += $6;
    nc++;
}
END{
    print xc/nc, yc/nc, zc/nc;
}

