#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
//#include <NTL/ZZ_pXModulus.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "CubicIdeal/cubic_ideal.h"
#include "Infideal/infrastructure_ideal.h"
#include "Hbounds/hbounds.h"
#include "Invariants/invariants.h"

using namespace std;

NTL_CLIENT

extern ZZ q;        // The characteristic.

extern ZZ_pX rho;            // rho = (x^2G(x))^(1/3)
extern ZZ_pX omega;          // omega = (xG^2(x))^(1/3)
extern int prec;
extern int belowinit;

extern ZZ_pX G;              
extern ZZ_pX H;
extern ZZ_pX f;  // f = GH^2. The function field is F_q(C), C: Y^3 = f.

extern int factordegs[21]; // The degrees of the prime factors of G and H.

int bsgs, kang;    // Run either Baby-Step Giant-Step or the Kangaroo Method.

extern int genus;     // The genus of the curve.
extern int rank;      // The unit rank.
extern int lambda;    // The bound used for computing E and L.
extern double alpha;  // An optimizer to adjust L2;
int irred;     // 1 if f is irreducible, 0 if not.
int manext;    // 0 to auto-extract R^S, 1 to manually extract R^S.  
long hashbits; // The number of bits in the size of the hashtable.

extern ZZ E;
extern ZZ U;
extern ZZ L;
ZZ order0; // A possible order.
ZZ order1; // A possible order.
ZZ N1 = to_ZZ(-1);

int maxCollision;
int type;


// Functions.

// Functions for Step 2: Determining extra information about the order h.

void smallOrder0(ZZ &a, ZZ &l);
int checkPrimes(ZZ &ord);
int shareFactors(ZZ order, ZZ factor);

// Functions for Step 3A: Baby-Step, Giant-Step

void bsgs0(ZZ& a, ZZ& factor);
void bsgs1();
void bsgs2();
//ZZ makeBabySteps(cubic_ideal &g, cubic_ideal &last, cubic_ideal* baby, ZZ* loc, long hashsize, ZZ& steps, ZZ& stepLen, int verbose);
ZZ makeBabySteps(cubic_ideal &g, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize, ZZ& steps, ZZ& a, ZZ& stepLen, int verbose);
ZZ makeBabySteps1(ZZVec& baby1, ZZVec& baby2, ZZVec& baby_deg, long hashsize, infrastructure_ideal& GS);
ZZ makeBabySteps2(ZZVec& baby1, ZZVec& baby2, ZZVec& baby_deg, long hashsize, infrastructure_ideal& GS);

// Functions for Step 3B: Pollard's Kangaroo

void kangaroo0(ZZ& a, ZZ& factor);
inline int vmap(cubic_ideal &A);
int distpoint(cubic_ideal &roo, int torw, ZZVec &dpt1, ZZVec &dpt2, ZZVec &dptd, ZZVec &dpw1, ZZVec &dpw2, ZZVec &dpwd, ZZ distance, ZZ &match, long size);

void kangaroo1();
inline int vmap(infrastructure_ideal &A);
int distpoint(infrastructure_ideal &roo, int torw, ZZVec &dpt1, ZZVec &dpt2, ZZVec &dptd, ZZVec &dpw1, ZZVec &dpw2, ZZVec &dpwd, ZZ &match, long size);
//int distinguished_point_check(infrastructure_ideal A);

// Functions for hashing and storage.
inline long roohash(cubic_ideal &A, long size);
inline void rooinsert(cubic_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, ZZ distance, long size);
inline void rooinsert(infrastructure_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, long size);
inline long roosearch(cubic_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, long hashsize);
inline long roosearch(infrastructure_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, long hashsize);
inline long hash(cubic_ideal &A, long size);
inline long hash(infrastructure_ideal &A, long size);
inline void insert(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, ZZ val, long size);
inline void insert(infrastructure_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long size);
inline long search(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize);
inline long search(infrastructure_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize);
inline ZZ baby1val(cubic_ideal &A);
inline ZZ baby2val(cubic_ideal &A);
inline ZZ baby1val(infrastructure_ideal &A);
inline ZZ baby2val(infrastructure_ideal &A);

// Functions for initialization.

int getinput(char *input);
void printInitData();
void printInitDataFile(char *file);
