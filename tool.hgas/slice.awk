#argv: mode (x:0 y:1 z:2), dx
BEGIN{
    dx = 1e30;
    dy = 1e30;
    dz = 1e30;
}
{
    if(ARGC < 3) {
        printf("awk -f slice.awk <dx=1e30> <dy=1e30> <dz=1e30> <datafile>\n") | "cat 1>&2";
        exit;
    }

    if(- dx < $4 && $4 < dx && - dy < $5 && $5 < dy && - dz < $6 && $6 < dz) {
        print $0;
    }
}