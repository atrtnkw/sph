{
    for(i = 0; i < 12; i++) {
	if(i == 0) {
	    printf("%2d %+.3e\n", 4, $8);
	} else {
	    printf("%2d %+.3e\n", 8+i*4, $(8+i));
	}
    }
}
