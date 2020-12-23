/*************************************************************************
*	lenstra.cc -- Use Lenstra's Elliptic Curves Algorithm to factor a 
*	large integer.  Works well when p+1+d is a product of small primes, 
*	where p is some factor of n, and d depends on the curve.  The algorithm 
*	was taken from "Rational Points on Elliptic Curves", by Silverman/Tate, 
*	pages 133-138.
*
*	Compile:
*	g++ -s -O4 -o lenstra lenstra.cc
*	Invoke:
*	./lenstra NumbertoFactor InitialK
*
* ChangeLog:
*  951124 -- Created by Ami Fischman <fischman@math.ucla.edu>
*  970204 to 970325 - minor bug fixes and cleaning up - Paul Herman <a540pau@pslc.ucla.edu>
**************************************************************************/

#include <Integer.h>

#define yes	1
#define no	0
#define K_incr 1	// K += K_incr  get new curve when old one bonks out on us
#define b_incr 30	// b += b_incr  when trying a new curve

#define max_tries_b  40	// number of b changes before changing K

// Initial point (x, y)  determines the shape of the curve
#define START_X  1	// x coordinate of the initail point to try
#define START_Y  1	// y coordinate of the initail point to try


/***************************
*   compute_k(Integer K)
*
*   choose k suitible for lenstra's alg.
*   Lots of people have different "good" k's.
****************************/
Integer compute_k(Integer K) {
	return LCM(K);
//	return fact(K);
}

/***********************
*  LenstraFactor(Integer n, Integer K)
*
*  find a factor using Lenstra's eliptical curve algorithm
************************/
Integer LenstraFactor(Integer n, Integer K) {
    Integer b, c, t, k, g, x, y, tempx, tempy, sum_x, sum_y, slope;
    int i, bits, newcurve, b_tries; 

    if (n == 6) return (Integer) 2;
    if (gcd(n,6) != 1) return gcd(n, 6);

    b = 1;					// setup our starting curve
    k = compute_k(K);
    bits = log2(k);				// i.e. log base 2
    b_tries = 0;
    for (;;) {
	x = START_X; y = START_Y;		// setup our starting points
        c = (y*y) - (x*x*x) - (b*x);		// calc. curve coefficients:


	// Make sure curve is OK to use (check Discriminant)
	t = gcd( 4*(b*b*b) + 27*(c*c), n);	// calc. Discriminant
	while (t != 1) {			// if 1<t<n, we got a factor!
	    if (t == n) {			// Discrim. must not == n
		b += b_incr;
		c = y*y - x*x*x - b*x;
		t = gcd( 4*(b*b*b) + 27*(c*c), n);
		continue;
	    }
	    return t;		// Hey, we got a factor along the way!
	}			// By now, gcd(Discr., n) = 1

	// Now we are ready to start adding points on the curve
	// This is done by calculating 2(x,y), 4(x,y), 8(x,y), ...
	// and adding them to (sum_x, sum_y) when needed.
	sum_x = 0; sum_y = 0; newcurve = no;

	if (b_tries > max_tries_b) {newcurve = yes;}
	else for (i=1; i<=bits; i++) {		// compute k*(x, y)
						// using powers of 2
	// Double the point {x, y}
	    tempx = x; tempy = 2*y;
	    g = gcd(tempy, n);
	    if (1<g && g<n) return g;
	    if (g == n) {newcurve = yes; break;}
	    slope = InvModN(tempy, n) * ( 3*(tempx*tempx) + b);	//now, g==1
	    x = (slope*slope - 2*tempx) % n;
	    y = (slope*(tempx-x) - y) % n;

	    if (!testbit(k, i)) continue;
	    if (!sum_y) {sum_x = x; sum_y = y; continue;}

	//Need to add {x,y} in to sum_{x,y}
	    tempx = (x - sum_x) % n;
	    tempy = (y - sum_y) % n;
	    g = gcd( tempx, n);
	    if (1<g && g<n) return g;
	    if (g == n) {newcurve = yes; break;}
	    slope = InvModN(tempx ,n) * tempy;		// at this point, g==1
	    tempx = ( (slope*slope - sum_x) - x)     % n;
	    tempy = ( slope*(sum_x - tempx) - sum_y) % n;
	    sum_x = tempx; sum_y = tempy;
	} /* for (i) */
	if (newcurve) {		// Give up and try a different curve & K
		b = 1;				// new constant
		b_tries = 0;			// reset iterations
		do {K += K_incr;} while (k == compute_k(K) );
		k = compute_k(K);		// new k
		bits = log2(k);
		}
	else {			// try a new b
		b += b_incr;	// new curve
		b_tries++;	// keep track how many times we do this
		}
    } /* while (1) */
}

Integer factor(Integer n, Integer K) {
	Integer f;
	f = LenstraFactor(n, K);
	if (isprime(f)) return f;
	else return factor(f, K);
}

Integer optimal_K(Integer n) {
	Integer K;
	int t;
	// try to guess optimal K, based on n

	K = logn(n, 10);
	K = K*K/2;
	if (K<2) K=2;
	return K;
}

int main(int argc, char *argv[]) {
    Integer n, K, f;
    int t;

    if (argc < 2) {
	cerr << "Usage: " << argv[0] << " NumberToFactor [InitialK]" << endl;
	cerr << "For an example, try " << argv[0] << " 1715761513 25" << endl;
	return -1;
    }

    n = argv[1];
    if (argc == 3) K = argv[2];
    else K = optimal_K(n);		// default, if none given

    cout << "Using initial K = " << K << endl;

    cout << n << ": " << flush;
    while (!isprime(n)) {
	f = factor(n, K);
	cout << f << " " << flush;
	n /= f;
	if (argc != 3) K = optimal_K(n);
	}
    cout << n << endl;
    return 0;
}
