# argv: xmax (ymax, zmax)
END{
    if(xmax == 0. && ymax == 0. && zmax == 0.) {
        print "Error: argv is required for xmax (ymax, zmax)" > "/dev/stderr"
        exit;
    }

    nmsh = 1000;
    rmax = sqrt(xmax**2 + ymax**2 + zmax**2);
    dmsh = rmax / nmsh;

    sim_gamma      = 5. / 3.;
    sim_rhoAmbient = 1.67e-24;
    sim_pAmbient   = 1.380658e-14;
    sim_vAmbient   = 0.;
    sim_expEnergy  = 1e51;
    sim_rInit      = 1e18;
    sim_vctr       = 4./3. * 3.141 * sim_rInit**3
    sim_mCenter    = 1.989e33;

    sim_rhoCenter  = sim_mCenter / sim_vctr;
    sim_uCenter    = sim_expEnergy / sim_mCenter;
    sim_pCenter    = sim_rhoCenter * (sim_gamma - 1.) * sim_uCenter;
    sim_vCenterCo  = sqrt(2. * sim_expEnergy / (3. * sim_mCenter)) / sim_rInit;

    for(i = 0; i < nmsh; i++) {
        r = i * dmsh;
        if(r < sim_rInit) {
            dens = sim_rhoCenter;
            pres = sim_pCenter;
            velr = sim_vCenterCo * r;
        } else {
            dens = sim_rhoAmbient;
            pres = sim_pAmbient;
            velr = 0.;
        }
#        printf("%+e %+e %+e %+e\n", r, dens, pres, velr);
         printf("%+e %+e %+e\n", dens, pres, velr);
    }
}

