# Not yet: mesh size
# argv: col=1,2
BEGIN {
    dx = 1e6;
}
{
    if(col == 1) {
        xcol = +2.225e8;
    } else if(col == 2) {
        xcol = -2.225e8;
    } else {
        printf("awk -f extractColumn1d.awk col=<1,2>\n") > "/dev/stderr";
        exit;
    }
    if(xcol-dx < $1 && $1 < xcol+dx) {
        print $0;
    }
}
