/***********************************************************************
*	shanks.cc - Uses shanks method to factor a large integer
*	Method found from some ubasic code out there somewhere.
*
*
*	Compile:
*	 g++ -s -O4 -o shanks shanks.cc
*	Invoke:
*	 ./shanks NumberToFactor
*	Where NumberToFactor is the number you wish to factor
*
* ChangeLog:
* 970307 -- Created by Paul Herman <a540pau@pslc.ucla.edu>
***********************************************************************/

#include <Integer.h>

Integer factor(Integer n)
{
    Integer a, f, h1, h2, k, p, pp, q, qq, qqq, r, te;
    int i = 0;

    k = sqrt(n);
    if (issquare(n)) return k;
    a=k; h1=k; h2=1; pp=0; qq=1; qqq=n; r=0;

    for (;;) {
	p = k-r;
	q = qqq + a*(pp-p);
	a = (p+k) / q; r = (p+k) % q;
	te = a*h1 + h2;
	h2 = h1;
	h1 = te;
	pp = p;
	qqq = qq; qq = q;
	te = sqrt(q);
	if ( (++i%2) != 0 || !issquare(q)) continue;
	te = h2-te;
	f = gcd(te, n);
	if (1<f && f<n) {
	    if (isprime(f)) return f;
	    else return factor(f);
	    }
	}
    return 0;
}


main (int argc, char *argv[])
{
	Integer n, t;
	int x0, a;
	int p;

	if (argc != 2) {
	   cerr << "Usage: " << argv[0] << " NumberToFactor" << endl;
	   return -1;
	   }

	n = argv[1];

	cout << n << ": " << flush;
	while (!isprime(n)) {
		t = factor(n);
		if (t==0) break;
		cout << t << " " << flush;
		n /= t;
		}
	cout << n << endl;
	return 0;
}
