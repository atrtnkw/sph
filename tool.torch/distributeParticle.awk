# argv: odir
{
    id=$2;
    out = sprintf("%s/id%010d.in", odir, id);
#    print $1, $3, $4 > out;
    print $1, $3, $4 >> out;
}
