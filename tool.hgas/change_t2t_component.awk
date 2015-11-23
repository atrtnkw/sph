# argv: He, C, O, Ne, Mg (fraction by mass)
{
    if(He + C + O + Ne + Mg != 1.) {
        print "The sum of He, C, O, Ne and Mg is not unity!" > "/dev/stderr";
        exit;
    }
    printf("%8d %2d %+e", $1, $2, $3);
    printf(" %+e %+e %+e", $4, $5, $6);
    printf(" %+e %+e %+e", $7, $8, $9);
    printf(" %+e %+e %+e %+e", $13, $14, $15, $17);
    printf(" %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e", He, C, O, Ne, Mg, $37, $38);
    printf(" %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e", $39, $40, $41, $42, $43, $44);
    printf("\n");
}
