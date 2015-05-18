###### argv: xmax, xmin, ymax, ymin, width
###### x: $3, y: $4
{
    if(NR == 1) {
        nx = int((xmax - xmin) / width);
        ny = int((ymax - ymin) / width);
        for(i = 0; i < nx; i++) {
            for(j = 0; j < ny; j++) 
                den[i,j] = 0.;
        }
    }

    if((xmin <= $3 && $3 < xmax) && (ymin <= $4 && $4 < ymax)) {
        ix = int(($3 - xmin) / width);
        iy = int(($4 - ymin) / width);
        den[ix,iy] += 1;
    }
}
END{
    denmax = 0.0;
    denmin = 1e30;
    for(i = 0; i < nx; i++) {
        for(j = 0; j < ny; j++) {
            if(den[i,j] > denmax)
                denmax = den[i,j];
            if(den[i,j] < denmin)
                denmin = den[i,j];
        }
    }

    denscale = 1. / (denmax - denmin);
    for(i = 0; i < nx; i++) {
        for(j = 0; j < ny; j++) {
            den[i,j] = (den[i,j] - denmin) * denscale;
        }
    }    

    for(i = 0; i < nx; i++) {
        for(j = 0; j < ny; j++) {
            x = i * width + xmin;
            y = j * width + ymin;
            printf("%+e %+e %+e\n", x, y, den[i,j]);
        }
        printf("\n");
    }
}
