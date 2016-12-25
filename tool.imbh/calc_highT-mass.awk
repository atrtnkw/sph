# argv: Tmaxcrit
BEGIN{
    m  = 0.;
    ms = 1.989e33;
}
{
    if($46 > Tmaxcrit) {
	m += $3;
    }
}
END{
    printf("%+e\n", m / ms);
}

