#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "cubic_ideal.h"
#include "Infideal/infrastructure_ideal.h"

using namespace std;

NTL_CLIENT

ZZ q;        // The characteristic.

extern ZZ_pX rho;            // rho = (x^2G(x))^(1/3)
extern ZZ_pX omega;          // omega = (xG^2(x))^(1/3)
extern int prec;

ZZ_pX D;       // D(x) = G(x)*H(x)^2
ZZ_pX G;              
ZZ_pX H;

int factordegs[9]; // The degrees of the prime factors of G and H.

int bsgs, kang;    // Run either Baby-Step Giant-Step or the Kangaroo Method.

int genus;     // The genus of the curve.
int rank;      // The unit rank.
int lambda;    // The bound used for computing E and L.
double alpha;  // An optimizer to adjust L2;
int irred;     // 1 if D is irreducible, 0 if not.
long hashbits; // The number of bits in the size of the hashtable.

long splitting[3];  // Stores the number of primes in each splitting category.
                    // 0: (P) = (p)^3 i.e. ramification.
                    // 1: (P) = (p)   i.e. inert
                    // 2: (P) = (p1)(p2)(p3) i.e. split completely
                    // (3:) (P) = (p1)(p2) i.e. partial splitting
                    // (4:) (P) = (p1)(p2)^2 i.e. wild ramification
                    // Note: If q = 1 (mod 3), there is no partial splitting
                    // Note: In char != 3, there is no wild ramification

vec_ZZ SvCache;  // Stores values of Sv(a) for 1 <= v <= lambda.
                 // If q = 1 mod 3, it stores Sv(1) and Sv(3).
                 // If q = 2 mod 3, it stores Sv(1) and Sv(2).


int MU[15] = {1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1};
int primes[45] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199};

RR Epr;
ZZ Ep;
ZZ E;
ZZ E1;
RR E2;
RR logE1;
RR logE2;
ZZ L;
ZZ L12;
ZZ L2;
ZZ L32;
ZZ order0; // A possible order.
ZZ order1; // A possible order.
ZZ N1 = to_ZZ(-1);

int maxCollision;
int type;


// Functions.

// Functions for Step 1: Computing E and L.

void approxh();
void Sv1(int v, int a);
void Sv2(int v, int a);
int mu(int a);
ZZ_p chi(ZZ_pX P);
ZZ_p rec3(ZZ_p a);
RR min(RR a, RR b);

// Functions for Step 2: Determining extra information about the order h.

ZZ smallOrder0();
int checkPrimes(ZZ &ord);
int shareFactors(ZZ order, ZZ factor);

// Functions for ranks 1 and 2.



// Functions for Step 3A: Baby-Step, Giant-Step

void bsgs0(ZZ& factor);
void bsgs1();
void bsgs2();
//ZZ makeBabySteps(cubic_ideal &g, cubic_ideal &last, cubic_ideal* baby, ZZ* loc, long hashsize, ZZ& steps, ZZ& stepLen, int verbose);
ZZ makeBabySteps(cubic_ideal &g, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize, ZZ& steps, ZZ& stepLen, int verbose);
ZZ makeBabySteps1(ZZVec& baby1, ZZVec& baby2, ZZVec& baby_deg, long hashsize, infrastructure_ideal& GS);

// Functions for Step 3B: Pollard's Kangaroo

void kangaroo0(ZZ& factor);
inline int vmap(cubic_ideal &A);
int distpoint(cubic_ideal &roo, int torw, ZZVec &dpt1, ZZVec &dpt2,ZZVec &dptd, ZZVec &dpw1, ZZVec &dpw2, ZZVec &dpwd, ZZ distance, ZZ &match, long size);


// Functions for hashing and storage.
inline long roohash(cubic_ideal &A, long size);
inline void rooinsert(cubic_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, ZZ distance, long size);
inline long roosearch(cubic_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, long hashsize);
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

int isCube(ZZ_p a);
int getinput(char *input);
void printInitData();
void printInitDataFile(char *file);

// Check the order.
int check(ZZ &order);
