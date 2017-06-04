if test $# -ne 2
then
    echo "sh $0 <ftype> <list>"
    exit
fi

ftyp=$1
list=$2

sort -n "$ftyp"_*.dat > tmp.snap

awk -f - flist="$list" tmp.snap<<EOF
{
    if(NR == 1) {
        n = 0;
        while((getline x < flist) > 0) {
            list[n] = int(x);
            n++;
        }
        list[n] = 1e9;
        i = 0;
    }
    if(\$1 > list[i]) {
        printf("# Lost ejected particle %8d!\n", list[i]);
        i++;
    } else if(\$1 == list[i]) {
        print \$0;
        i++;
    }
}
EOF

rm -f tmp.*

