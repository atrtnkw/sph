BEGIN{
    eb = 0.;
    ek = 0.;
    ep = 0.;
    n  = 0;
}
{
    eb+=0.5*$3*($7**2+$8**2+$9**2+$25);
    ek+=0.5*$3*($7**2+$8**2+$9**2);
    ep+=0.5*$3*$25;
    n++;
}
END{
    printf("n: %8d eb: %+e ek: %+e ep: %+e\n", n, eb, ek, ep)
}
