BEGIN{
    n = 0;
    while((getline var < "ref/triplealpha.dat") > 0) {
        dens[n] = substr(var,  1, 12) + 0.;
        temp[n] = substr(var, 14, 12) + 0.;
        n++;
    }
    dens[n] = 1.582045e+05;
    temp[n] = 1e10;
}
{
    dens0 = $16;
    temp0 = $21;
    for(i = 0; i < n; i++) {
        if(temp0 < temp[i]) {
            break;
        } else if(temp[i] <= temp0 && temp0 < temp[i+1]) {
            if(dens0 >= dens[i+1]) {
                print $0;
            }
            break;
        }
    }
}

