# argv: col=1,2,-1, dx, px
BEGIN {
    xcol = 1e30;
}
{
    if(col == 1) {
        if(xcol != 1e30) {
            printf("awk -f extractColumn1d.awk col=1 dx=???\n") > "/dev/stderr";
            exit;
        }
        xcol = +2.225e8;
    } else if(col == 2) {
        if(xcol != 1e30) {
            printf("awk -f extractColumn1d.awk col=2 dx=???\n") > "/dev/stderr";
            exit;
        }
        xcol = -2.225e8;
    } else if (col == -1) {
        if(xcol == 1e30) {
            printf("awk -f extractColumn1d.awk col=-1 xcol=??? dx=???\n") > "/dev/stderr";
            exit;
        }
    } else {
        printf("awk -f extractColumn1d.awk col=<-1,1,2> dx=???\n") > "/dev/stderr";
        exit;
    }
    if(dx == 0) {
        printf("awk -f extractColumn1d.awk col=<-1,1,2> dx=???\n") > "/dev/stderr";
        exit;
    }
    if(xcol-dx < $1 && $1 < xcol+dx) {
        print $0;
    }
}
