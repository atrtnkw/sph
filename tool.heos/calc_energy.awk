BEGIN{
    etot = 0.0;
    ekin = 0.0;
    epot = 0.0;
    ethm = 0.0;
}
{
    mass = $2;
    xvel = $6;
    yvel = $7;
    zvel = $8;
    uene = $12;
    pot  = $24;

    ikin = 0.5 * mass * (xvel**2 + yvel**2 + zvel**2);
    ipot = 0.5 * mass * pot;
    ithm =       mass * uene;

    ekin += ikin;
    epot += ipot;
    ethm += ithm;
    etot += (ikin + ipot + ithm);
}
END{
    printf("%+e %+e %+e %+e\n", etot, ekin, epot, ethm);
}
