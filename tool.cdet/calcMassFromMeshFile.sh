if test $# -ne 2
then
    echo "sh $0 <ifile> <iline>" >&2
    echo "iline= 8: 4He"  >&2
    echo "iline= 9: 12C"  >&2
    echo "iline=10: 16O"  >&2
    echo "iline=11: 20Ne" >&2
    echo "iline=12: 24Mg" >&2
    echo "iline=13: 28Si" >&2
    echo "iline=14: 32S"  >&2
    echo "iline=15: 36Ar" >&2
    echo "iline=16: 40Ca" >&2
    echo "iline=17: 44Ti" >&2
    echo "iline=18: 48Cr" >&2
    echo "iline=19: 52Fe" >&2
    echo "iline=20: 56Ni" >&2
    exit
fi

ifile=$1
iline=$2

awk 'BEGIN{pr=0}{print 4/3.*3.14*($1**3-pr**3)*$2*$(iline);pr=$1}' iline=$iline $ifile \
    | awk 'BEGIN{sum=0}{sum+=$1}END{printf("%+e\n", sum/1.989e33)}'
