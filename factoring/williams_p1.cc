/**************************************************************************
*	williams_p1.cc - Uses William's p+1 method to factor large integer
*	Method found from some ubasic code out there somewhere.
*
*
*	Compile:
*	g++ -s -O4 -o shanks shanks.cc
*	Invoke:
*	./shanks NumberToFactor
*	Where NumberToFactor is the number you wish to factor
*
* ChangeLog
*  970307 -- Created by Paul Herman <a540pau@pslc.ucla.edu>
**************************************************************************/


#include <Integer.h>

Integer factor(Integer n, Integer p)
{
    Integer d, count, vs, vb;
    int bp, cc, i, k;

    cc = 10; count = 0;

    for (;;) {
	d = gcd(p-2, n);
	if (1<d && d<n) {
	    if (isprime(d)) return d;
	    else factor(d, p);
	    }
	for (i=1; i<=cc; i++) {
	    count++;
	    vs = p;
	    vb = (p*p - 2)%n;
	    bp = log2(count);
	    for (k=bp; k>-1; k--) {
		if (testbit(count, k)) {vs = (vs*vb - p)%n; vb = (vb*vb - 2)%n;}
		else {vb = (vs*vb - p)%n; vs = (vs*vs - 2)%n;}
		}
	    p = vs;
	    }
	}
    return 0;  // pray this never happens...
}


main (int argc, char *argv[])
{
	Integer n, t, p;
	int x0, a;

	if (argc < 2) {
	   cerr << "Usage: " << argv[0] << " NumberToFactor [p value]" << endl;
	   return -1;
	   }

	n = argv[1];
	p = 4;
	if (argc == 3) p = argv[2];

	cout << n << ": " << flush;
	while (!isprime(n)) {
		t = factor(n, p);
		if (t==0) break;
		cout << t << " " << flush;
		n /= t;
		}
	cout << n << endl;
	return 0;
}
