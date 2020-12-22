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

using namespace std;

NTL_CLIENT

ZZ q;        // The characteristic.

int degG, degH;  // The degrees of G and H.

ZZ_pX D;       // D(x) = G(x)*H(x)^2
ZZ_pX G;              
ZZ_pX H;

int factordegs[20]; // The degrees of the prime factors of G and H.

int genus;     // The genus of the curve.
int rank;      // The unit rank.
int lambda;    // The bound used for computing E and L.
double alpha1, alpha2, alpha3;  // An optimizer to adjust L2;
int irred;     // 1 if D is irreducible, 0 if not.
long hashbits; // The number of bits in the size of the hashtable.

vec_ZZ SvCache;  // Stores values of Sv(a) for 1 <= v <= lambda.
                 // If q = 1 mod 3, it stores Sv(1) and Sv(3).
                 // If q = 2 mod 3, it stores Sv(1) and Sv(2).

time_t set1, time1;
time_t setb[3], timeb[3];
time_t setg[3], timeg[3];


int MU[15] = {1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1};
int primes[45] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199};

RR psi1, psi2, psi3;
ZZ E;
ZZ Es[3];
RR E1r, E2r;
RR logE1;
RR logE2;
ZZ L, Ls[3];
ZZ L2[3];
int best;
ZZ h;
ZZ HasseLow;
ZZ HasseHigh;
ZZ order0; // A possible order.
ZZ order1; // A possible order.
ZZ N1 = to_ZZ(-1);

int samples;
int maxCollision;

ZZ babysteps[3];
ZZ giantsteps[3];

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
int checkPrimes(ZZ& ord);
int shareFactors(ZZ order, ZZ factor);

// Functions for Step 3: Baby-Step, Giant-Step

void bsgs0(ZZ& factor, int run);
void bsgs1();
void bsgs2();

ZZ makeBabySteps(cubic_ideal &g, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize, ZZ& steps, ZZ& stepLen, int verbose);
inline long hash(cubic_ideal &A, long size);
inline void insert(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, ZZ val, long size);
inline long search(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize);
inline ZZ baby1val(cubic_ideal &A);
inline ZZ baby2val(cubic_ideal &A);

// Functions for initialization.

int isCube(ZZ_p a);
int getinput(char *input);
void setConstants();
void getPolynomial();
