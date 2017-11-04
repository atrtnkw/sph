# flist: txxxx_x.log.id
{
    if(NR == 1) {
        n = 0;
        while((getline x < flist) > 0) {
            list[n] = int(x);
            n++;
        }
#        for(i = 0; i < n; i++) {
#            print list[i];
#        }
        ifirst = 0;
    }    
    for(i = ifirst; i < n; i++) {
        if($1 == list[i]) {
            print $1;
            ifirst = i;
        }
    }
    if(NR % 10000 == 0) {
        print $1 > "/dev/stderr";
    }
}
