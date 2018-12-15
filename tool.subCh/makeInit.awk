BEGIN{
    # base parameter
    dnbs = 1e7;
    tpbs = 1e8;
    vxbs = 0.;
    # hotspot parameter
    xdis = 1e7;
    tdis = 3e9;
    
    xmax = 1e8;
    nmax = 800;

    dx   = 2 * xmax / nmax;
    for(i = 0; i < nmax; i++) {
        px = dx * i - xmax;
        if(px*px < xdis*xdis) {
            dn = dnbs;
            vx = vxbs;
            xx = (px > 0.) ? px : -px;
            tp = tdis - (tdis - tpbs) * (xx / xdis) + tpbs;
        } else {
            dn = dnbs;
            vx = vxbs;
            tp = tpbs;
        }
        printf("%+e %+e %+e %+e\n", px, dn, tp, vx);
    }
}
