#if test $# -ne 2
#then
#    echo "sh $0 <idir0> <idir1>"
#    exit
#fi

echo "1" > hoge.dum

for ndir in $(seq -f "%03g" 191 1 218)
do
    idir=t"$ndir"
#    idir=dummy"$ndir"
    flag=0
    for time in $(seq -f "%06g" 0 1 999999)
    do
        ifile="$idir"/t"$time"_p000000.hexa
        if test $flag -eq 0
        then
            if test -e  $ifile
            then
                flag=1
            del0=`awk '{printf("%d\n", time+1);}' time=$time hoge.dum`
            fi
        else
            if ! test -e $ifile
            then
                del1=`awk '{printf("%d\n", time-2);}' time=$time hoge.dum`
                break
            fi
        fi
    done
    
    for time in $(seq -f "%06g" $del0 1 $del1)
    do
        rm -f "$idir"/t"$time"_p*.hexa
    done
done
    
rm -f hoge.dum
