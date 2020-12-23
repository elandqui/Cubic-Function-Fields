#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <ctime>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "../CubicIdeal/cubic_ideal.h"
#include "../Infideal/infrastructure_ideal.h"
#include "../Hbounds/hbounds.h"
#include "../Invariants/invariants.h"

using namespace std;

NTL_CLIENT

extern ZZ q;        // The characteristic.


extern ZZ_pX f;       // f(x) = G(x)*H(x)^2
extern ZZ_pX G;              
extern ZZ_pX H;

extern ZZ One;
extern ZZ Zero;

extern int factordegs[21]; // The degrees of the prime factors of G and H.

extern int genus;     // The genus of the curve.
extern int rank;      // The unit rank.
extern int lambda;    // The bound used for computing E and U.
extern double alpha;  // An optimizer to adjust L;
double tauratio;              // A rounded GS/BS ratio for optimization.
int m;         // The number of kangaroos.
int mt, mw;    // The number of tame and wild kangaroos.
int irred;     // 1 if D is irreducible, 0 if not.
int phase1time;

extern long splitting[4];  // Stores the number of primes in each splitting category.
                    // 0: (P) = (p)^3 i.e. ramification.
                    // 1: (P) = (p)   i.e. inert
                    // 2: (P) = (p1)(p2)(p3) i.e. split completely
                    // 3: (P) = (p1)(p2) i.e. partial splitting
                    // (4:) (P) = (p1)(p2)^2 i.e. wild ramification
                    // Note: If q = 1 (mod 3), there is no partial splitting
                    // Note: In char != 3, there is no wild ramification

extern vec_ZZ SvCache;  // Stores values of Sv(a) for 1 <= v <= lambda.
                 // If q = 1 mod 3, it stores Sv(1) and Sv(3).
                 // If q = 2 mod 3, it stores Sv(1) and Sv(2).

int distbits;    // The number of zero bits to determine 
                       // if an ideal is distinguished.
vec_ZZX dp;    // Storage for distinguished points.
vec_ZZ dists;    // The distances.
char *type;
int *roonum;
int size;        // Size of the distinguished point arrays.

extern int MU[15];// = {1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1};
extern int primes[45];// = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199};

extern ZZ E;
extern ZZ L;
extern ZZ U;

ZZ order;
ZZ R;

/**************/
/* Functions. */
/**************/

// Functions for Step 1: Computing E and L.

void approxh();

// Functions for Step 3B: Pollard's Kangaroo and processing.
void getRooJumps0(ZZ &a, ZZ &l);
void getRooJumps1();
void getRooJumps2(ZZ &a, ZZ &l);
int readDPfile(ZZ &l);
void sort(int hi, int lo);
int search(int len, int &matchloc, ZZ &l);
void reassign(char tw, int num, ZZ &l);

// Extra functions for the unit rank 2 case.
//void latticeRoo(ZZ &icn, ZZ &R, ZZ &e10, ZZ &e12, ZZ &e20, ZZ &e22);

// Functions for initialization.
int getinput(char *input);
