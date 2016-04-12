# argv: mode (0: bound, 1: unbound)
{
    eb = 0.5*($7**2+$8**2+$9**2)+$45;
    if(mode == 0) {
        if(eb < 0.) {
            print $0;
        }
    } else {
        if(eb >= 0.) {
            print $0;
        }
    }
}
