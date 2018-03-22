# reduction rate: n (=2^integer)
{
    if(NR == 1) {
        print $0;
    } else {
        if(NR % n == 1) {
            $3  *= n;
            $13 *= n;
            print $0;
        }
    }
}
