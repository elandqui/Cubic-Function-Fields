#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pXFactoring.h>
#include <iostream>
#include <NTL/RR.h>

using namespace std;

NTL_CLIENT

extern ZZ q;        // The characteristic.

extern ZZ_pX G;     // The defining polynomials.
extern ZZ_pX H;     // The defining polynomials.
extern ZZ_pX f;     // f = GH^2. The function field is F_q(C), C: Y^3 = f.
extern int genus;   // The genus of the curve.
extern int rank;    // The unit rank.

// Functions for Step 1: Computing E, U, and L.

void approxh(int verbose);
void Sv1(int v, int a, ZZ &start, ZZ &end);
void Sv2(int v, int a);
int mu(int a);
ZZ_p chi(const ZZ_pX &P);
ZZ_pX chi2(const ZZ_pX &P, const ZZ &Q);
ZZ_p rec3(ZZ_p a);
ZZ_p rec32(ZZ_p a, ZZ &Q);

// Functions for Step 2: Determining extra information about the order h.

void smallOrder0(ZZ &a, ZZ &l);
int checkPrimes(ZZ &ord);
int shareFactors(ZZ order, ZZ factor);

// Functions for initialization.

int isCube(ZZ_p a);
