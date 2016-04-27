#argv: px, py, pz, dx, dy, dz
BEGIN{
    px = 0.;
    py = 0.;
    pz = 0.;
    dx = 1e30;
    dy = 1e30;
    dz = 1e30;
}
{
    if(ARGC < 3) {
        printf("awk -f slice.awk <dx=1e30> <dy=1e30> <dz=1e30> <datafile>\n") | "cat 1>&2";
        exit;
    }

    if(px - dx < $4 && $4 < px + dx && py - dy < $5 && $5 < py + dy && pz - dz < $6 && $6 < pz + dz) {
        print $0;
    }
}