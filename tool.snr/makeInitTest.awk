# argv: xmax (ymax, zmax)
END{
    if(xmax == 0. && ymax == 0. && zmax == 0.) {
        print "Error: argv is required for xmax (ymax, zmax)" > "/dev/stderr"
        exit;
    }

    nmsh = 1000;
#    xmax = 2e19;
#    ymax = 0.;
#    zmax = 0.;
    rmax = sqrt(xmax**2 + ymax**2 + zmax**2);
    dmsh = rmax / nmsh;

    sim_gamma      = 5. / 3.;
    sim_rhoAmbient = 1.67e-24;
    sim_pAmbient   = 1.380658e-14;
    sim_vAmbient   = 0.;
    sim_expEnergy  = 1e51;
    sim_rInit      = 1e18;
    sim_vctr       = 4./3. * 3.141 * sim_rInit**3

    for(i = 0; i < nmsh; i++) {
        r = i * dmsh;
        if(r < sim_rInit) {
            sim_rhoCenter = sim_rhoAmbient;
            sim_pCenter   = (sim_gamma-1.) * sim_expEnergy / sim_vctr;
            sim_vCenter   = sim_vAmbient
#            printf("%+e\n", r);
            printf("%+e %+e %+e\n", sim_rhoCenter, sim_pCenter, sim_vCenter);
        } else {
#            printf("%+e\n", r);
            printf("%+e %+e %+e\n", sim_rhoAmbient, sim_pAmbient, sim_vAmbient);
        }
    }
}

