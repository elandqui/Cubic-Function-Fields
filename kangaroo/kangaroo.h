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

extern ZZ_pX f;
extern ZZ_pX G;              
extern ZZ_pX H;

extern int genus;
extern int rank;

extern ZZ One;
extern ZZ Zero;

//int tauratio;

int m;              // The number of kangaroos.
extern int lambda;  // The polynomial degree to examine.
extern long splitting[4];  // Stores the number of primes in each splitting category.
                    // 0: (P) = (p)^3 i.e. ramification.
                    // 1: (P) = (p)   i.e. inert
                    // 2: (P) = (p1)(p2)(p3) i.e. split completely
                    // 3: (P) = (p1)(p2) i.e. partial splitting
extern vec_ZZ SvCache;  // Stores values of Sv(a) for 1 <= v <= lambda.
                 // If q = 1 mod 3, it stores Sv(1) and Sv(3).
                 // If q = 2 mod 3, it stores Sv(1) and Sv(2).

char rootype;   // t if tame; w if wild.
int roonum;     // The number of the kangaroo (#roonum of m/2). 0 <= roonum < m/2
char *lastdpname; // The name of the file storing the last trap for this roo.
char *alldpname; // The name of the file storing all traps for this roo.

extern int factordegs[21]; // The degrees of the prime factors of G and H.

ZZ start;     // Starting position.
ZZ offset;    // Space the kangaroos offset units apart initially
int distbits; // Number of 0's to make a distinguished point.

// Kangaroo jumps and jump distances.
ZZ jumpdistance[64];
ZZ jumpdistance1[64];
ZZX gz[6];
ZZX jumpsz[384];

// Functions.
void getSplitting();
void kangaroo0();
void kangaroo1();
void kangaroo2();
void kangaroo2(ZZ &a, ZZ &l);
inline int vmap(cubic_ideal &A);
inline int vmap(infrastructure_ideal &A);
int distpoint(cubic_ideal &roo, ZZ &distance, long &dps, ZZ &jumps);
int distpoint(infrastructure_ideal &roo, long &dps, ZZ &bsjumps, ZZ &gsjumps, ZZ &eclass, ZZ &l);
void getRooInfo();

