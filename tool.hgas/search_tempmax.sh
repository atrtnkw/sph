if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>"
    exit
fi

idir=$1
tbgn=$2
tend=$3

for time in $(seq $tbgn 1 $tend)
do
    ptim=`printf "%04d" $time`
    ifile="$idir"/sph_t"$ptim".dat
    printf "%5d" $time
    awk 'BEGIN{rhol=10**5.0;rhoh=10**5.5;}{if(rhol<=$16&&$16<rhoh) print $21;}' $ifile \
        | awk 'BEGIN{tmax=0.;}{if($1>tmax){tmax=$1;}}END{printf(" %+e", tmax);}' t=$time
    awk 'BEGIN{rhol=10**5.5;rhoh=10**6.0;}{if(rhol<=$16&&$16<rhoh) print $21;}' $ifile \
        | awk 'BEGIN{tmax=0.;}{if($1>tmax){tmax=$1;}}END{printf(" %+e", tmax);}' t=$time
    awk 'BEGIN{rhol=10**6.0;rhoh=10**6.5;}{if(rhol<=$16&&$16<rhoh) print $21;}' $ifile \
        | awk 'BEGIN{tmax=0.;}{if($1>tmax){tmax=$1;}}END{printf(" %+e", tmax);}' t=$time
    awk 'BEGIN{rhol=10**6.5;rhoh=10**7.0;}{if(rhol<=$16&&$16<rhoh) print $21;}' $ifile \
        | awk 'BEGIN{tmax=0.;}{if($1>tmax){tmax=$1;}}END{printf(" %+e", tmax);}' t=$time
    awk 'BEGIN{rhol=10**7.0;rhoh=10**7.5;}{if(rhol<=$16&&$16<rhoh) print $21;}' $ifile \
        | awk 'BEGIN{tmax=0.;}{if($1>tmax){tmax=$1;}}END{printf(" %+e", tmax);}' t=$time
    printf "\n"
done
