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
#include <mpi.h>
#include "../CubicIdeal/cubic_ideal.h"
#include "../Infideal/infrastructure_ideal.h"
#include "../Hbounds/hbounds.h"
#include "../Invariants/invariants.h"

using namespace std;

NTL_CLIENT

extern ZZ q;           // The characteristic.

extern ZZ_pX f;       // f(x) = G(x)*H(x)^2
extern ZZ_pX G;              
extern ZZ_pX H;

extern ZZ One;
extern ZZ Zero;

int irred;      // 1 if D is irreducible, 0 if not.

int m;          // The number of kangaroos.
int mt, mw;     // The number of tame and wild kangaroos.

extern int genus;      // The genus of the curve.
extern int rank;       // The unit rank.
extern int lambda;     // The polynomial degree to examine.
extern double alpha;   // An optimizer to adjust L2;

ZZ a, l;        // If it is known that h = a mod l.
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

extern ZZ E;            // The approximation of the order.
extern ZZ L;
extern ZZ U;          // The error bound of the approximation.

ZZ order;        /* The order */
vec_ZZ factors;
vec_long exponents;

vec_ZZX dp;      // Storage for distinguished points.
vec_ZZ dists;    // The distances.
char *type;
int *nums;
int size;        // Size of the distinguished point arrays.

char rootype;    // t if tame; w if wild.
int roonum;      // The number of the kangaroo (#roonum of m/2). 0 <= roonum < m/2
char *lastdpname;   // The name of the distinguished point file for this roo.
char *alldpname;  // Where the roo stores all DPS.

int factordegs[21]; // The degrees of the prime factors of G and H.

extern int MU[15];
extern int primes[45];

ZZ start;     // Starting position.
ZZ offset;    // Space the kangaroos offset units apart initially
int distbits; // Number of 0's to make a distinguished point.

// Kangaroo jumps and jump distances.
ZZ jumpdistance[64];
ZZ jumpdistance1[64];
ZZX gz[6];
ZZX jumpsz[384];

/*************/
/* Functions */
/*************/

/********************************************/
/* Functions for Step 1: Computing E and L. */
/********************************************/

// For the master:

void approxh();
void Sv1(int v, int a);
void Sv2(int v, int a);

// For the kangaroos:

void getSplitting();

/****************************************************/
/* Functions for Step 2: Extra information about h. */
/****************************************************/


/************************************************************/
/* Functions for Step 3: Pollard's Kangaroo and processing. */
/************************************************************/

// For the master:

void getRooJumps0();
void getRooJumps1();
void getRooJumps2();
int readDPfile();
void sort(int hi, int lo);
int search(int len, int &matchloc);
void reassign(char tw, int num);

// For the kangaroos:

void kangaroo0();
void kangaroo1();
void kangaroo2();
inline int vmap(cubic_ideal &A);
inline int vmap(infrastructure_ideal &A);
int distpoint(cubic_ideal &roo, ZZ &distance, long dps, ZZ &jumps);
int distpoint(infrastructure_ideal &roo, long dps, ZZ &bsjumps, ZZ &gsjumps);

/********************************/
/* Functions for initialization */
/********************************/

int getRooInfo();
int getinput();
