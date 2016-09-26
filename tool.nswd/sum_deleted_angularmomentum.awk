#argv: tp
BEGIN{
    RS = "";
}
{
    if($1 != "###") {
	next;
    }

    t  = $3;

    if(t <= tp) {
	next;
    }

    id = $5;
    ms = $7;
    px = $56;
    py = $57;
    pz = $58;
    vx = $60;
    vy = $61;
    vz = $62;

    dangx = ms * (py * vz - pz * vy);
    dangy = ms * (pz * vx - px * vz);
    dangz = ms * (px * vy - py * vx);

    printf("%18.10f %8d %+.10e %+.10e %+.10e\n", t, id, dangx, dangy, dangz);
    
}

