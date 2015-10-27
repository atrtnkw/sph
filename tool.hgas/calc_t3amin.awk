BEGIN{
#    if($2 == 1) {
#        dens = $16;
#        temp = $21;
    gpm  = 4.; # helium [g/mol]
    cp   = 2.1e8/gpm;
    dens = 1e6;
    temp = 3e8;
#    e3a  = 5.09e11 * dens**2 / (temp*1e-8)**3 * exp(2.76e-3*dens**0.5/(temp*1e-8)**1.5 - 44.027/(temp*1e-8));
    e3a  = 5.09e11 * dens**2 / (temp*1e-8)**3 * exp(2.76e-3*dens**0.5/(temp*1e-8)**1.5 - 44.027/(temp*1e-8));
    t3a  = cp * temp / e3a;
    tdy  = 1. / (24 * 3.141 * 6.67259e-8 * dens)**0.5;
    print t3a, tdy;
#}
}
