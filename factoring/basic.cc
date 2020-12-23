/*************************************************
*	fact.cc -- Factor any integer into primes
*
*	This code actualy contains *two* factorization routines, both based on 
*	the same "try each number before sqrt(n)" approach.  They are:
*
*	smart_factor()
*	Like factor() except it skips all multiples of 2, 3 & 5 cutting the time 
*	in more than 1/2.  Just a tinge more "practical" than.... 
*
*	factor()
*	Uses simplistic "up to sqrt(n)" approach.  Not "practical".
*
*	Compile: 
*	g++ -s -O4 -o fact fact.cc
*	Invoke:
* 	./fact NumberToFactor
*
*  ChangeLog
*  950516 -- Created by Ami Fischman <fischman@math.ucla.edu>
*  970307 -- Added smart_factor() -- Paul Herman <a540pau@pslc.ucla.edu>
********************************************************/

#include <Integer.h>

int add[] = {4, 2, 4, 2, 4, 6, 2, 6};  // start after the number 7
// The next step: From 7, 7+add1[0] = 7+4 = 11, 11 + add1[1] = 11+2 = 13.
// Once through yields all primes up to 217 with no composites.
// Next iteration with add2 yields all primes up to 427, 
// but with some composites, etc.
// This skips all composites that are multiples of 2, 3, 5, 7, and 11.
// Length of add1[] = 44, add2[] = 45, 
// int add1[] = {4, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 14, 4, 6, 2, 10, 2, 6, 6, 4, 6, 6, 2, 10, 2, 4, 2, 12, 6};
// int add2[] = {4, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 14, 4, 6, 2, 10, 2, 6, 6, 4, 2, 4, 6, 2, 10, 2, 4, 2, 12, 6};

/***************************
*	smart_factor(Integer n, start)
*	Given n, try to succesivly find a factor.  Begin with start.
*	This is called "smart" because it skips multiples of 2, 3, & 5
****************************/
Integer smart_factor(Integer n, Integer start) {
	Integer f_try;
	int ai;

	if (n%2 == 0) return 2;
	if (n%3 == 0) return 3;
	if (n%5 == 0) return 5;
	if (n%7 == 0) return 7;

	f_try = start;
	if (f_try<7) f_try = 7;
	ai = Int2i((f_try-7)%8);	// align ai to the correct index

	for (;;) {
	   if (n%f_try == 0) {
		if (isprime(f_try)) return f_try;
		else return smart_factor(f_try, 7);
		}
	   else {
		f_try += add[ai++];
		ai %= 8;
		}
	   }
}

/***************************
*	factor(Integer n, start)
*	Given n, try to succesivly find a factor.  Begin with start.
****************************/
Integer factor(Integer n, Integer start) {
    Integer i;

    i = start;
    while (i < n) {
	if ( n%i == 0 ) {
	    if (isprime(i)) return i;
	    else return factor(i, 2);	// start could have been too high
	    }
	i++;
	}
    return 1;	// start could have been too high
};

int main(int argc, char *argv[]) {

    Integer n, f;
    char ai, *p, smart;

    if (argc < 2) {
	cerr << "Usage: " << argv[0] << " [-s] NumberToFactor" << endl;
	return 1;
    }

    // check for options
    ai = 1; smart = 0;
    p = (char *) argv[1];
    if (*p == '-') {
	switch (*(p+1)) {
	    case 's':
		smart = 1;
		break;
	    default:
		break;
	    }
	ai++;
	}
    
    n = argv[ai];
    cout << n << ": " << flush;

    if (isprime(n)) {cout << "is prime..." << endl; return 0;}

    f = 2;
    while (!isprime(n)) {
	if (smart) f = smart_factor(n, 7);
	else       f = factor(n, 2);
	cout << f << " " << flush;
	n /= f;
	}
    cout << n << endl;
    return 0;
};
