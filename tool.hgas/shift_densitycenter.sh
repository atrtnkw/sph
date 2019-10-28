if test $# -ne 1
then
    echo "sh $0 <file>" >&2
    exit
fi

file=$1
colden=16

awk 'BEGIN{xc=yc=zc=0.;dn=0.;}
{xc+=$(colden)*$4;yc+=$(colden)*$5;zc+=$(colden)*$6;dn+=$(colden)}
END{printf("%+.16e %+.16e %+.16e\n", xc/dn, yc/dn, zc/dn);}' colden=$colden $file > tmp.xc

xc=`awk '{print $1}' tmp.xc`
yc=`awk '{print $2}' tmp.xc`
zc=`awk '{print $3}' tmp.xc`

awk '{$4-=xc;$5-=yc;$6-=zc; print $0;}' xc=$xc yc=$yc zc=$zc $file
