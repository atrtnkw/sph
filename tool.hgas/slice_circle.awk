#argv: R, dR
BEGIN{
    R  = 0.;
    dR = 1e30;
}
{
    if(NR == 1) {
        printf("# R=%e dR=%e\n", R, dR);
        Rl2 = (R - dR)**2;
        Ru2 = (R + dR)**2;
    }

    Ri2 = $4**2 + $5**2;

    if(Rl2 < Ri2 && Ri2 < Ru2) {
        theta = atan2($5, $4);
        Rth   = R * theta;        
        print $0, Rth;
    }
}
