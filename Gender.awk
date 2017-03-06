BEGIN {
    if (x+y < 1000){
        print "UNKNOWN"
    } else if ((y/ytotal) / (x/xtotal) < 0.05) {
        print "FEMALE"
    } else if ((y/ytotal) / (x/xtotal) > 0.1) {
        print "MALE"
    } else {
	    print "UNKNOWN"
    }
}