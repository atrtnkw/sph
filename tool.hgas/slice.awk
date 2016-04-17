#argv: mode (x:0 y:1 z:2), dx
{
    if(ARGC != 4) {
        printf("awk -f slice.awk <mode=0/1/2> <dx=> <datafile>\n") | "cat 1>&2";
        exit;
    }

    if(mode == 0) {
        col = 4;
    } else if(mode == 1) {
        col = 5;
    } else {
        col = 6;
    }
    
    if(- dx < $(col) && $(col) < dx) {
        print $0;
    }
}