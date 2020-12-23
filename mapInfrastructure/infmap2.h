#include <NTL/ZZ_pX.h>
#include <iostream>
#include <time.h>
using namespace std;

#define MAX(a, b) (a > b ? a : b)

NTL_CLIENT

long P;      // The characteristic.
int u;      // A primitive cube root of 1 in F_p.

// We work in the function field F_p(x)[y]/<y^3 - x^2G(x)>, 
// where deg(G(x)) = 4. 
// The basis elements rho and omega are elements of F_p((1/x)),
// so they must be truncated and shifted to fit as a polynomial.

int shift=5;

ZZ_pX rho;            // rho = (x^2G(x))^(1/3)
ZZ_pX omega;          // omega = (xG^2(x))^(1/3)

ZZ_pX mu0, mu1, mu2;  // mu = mu0+mu1*rho+mu2*omega
ZZ_pX nu0, nu1, nu2;  // nu = nu0+nu1*rho+nu2*omega
ZZ_pX t0, t1, t2;     // theta = t0+t1*rho+t2*omega
ZZ_pX sig0, sig1, sig2;
ZZ_pX tau0, tau1, tau2;
ZZ_pX phi0, phi1, phi2;


ZZ_pX d;              // d = denominator of mu and nu
int ddeg=0;           // degree of d
ZZ_pX G;              // as above
ZZ_pX F;              // F(x) = G(x)*x^2

vec_ZZ_pX I;          // C stores theta and I stores mu and nu.
vec_ZZ_pX D;          // D stores the denominators.

int p, l;             // p = preperiod length, l = period length 
int m;                // m = period length of second cycle
int num;              // Number of infrastructure elements discovered.

int coordl[180];      // Stores coordinates of the infrastructure.
int coordm[180];
int deg0[180];        // Stores the degrees of the infrastructure.
int deg1[180];
int deg2[180];
int d0=0, d1=0, d2=0; // Current degrees of the infrastructure.
int dege1=0;          // Degrees of the units.
int dege2=0;          

ZZ_pX cuberoot(ZZ_pX A, int prec);
void basis(int prec);
inline int zetadeg(ZZ_pX a, ZZ_pX b, ZZ_pX c);
inline int xideg(ZZ_pX b, ZZ_pX c);
inline int etadeg(ZZ_pX b, ZZ_pX c);
inline int xideg1(ZZ_pX b, ZZ_pX c);
inline int etadeg1(ZZ_pX b, ZZ_pX c);
inline int xideg2(ZZ_pX b, ZZ_pX c);
inline int etadeg2(ZZ_pX b, ZZ_pX c);
inline int discdeg(ZZ_pX mu1,ZZ_pX mu2, ZZ_pX nu1, ZZ_pX nu2);
inline int eltdeg(ZZ_pX a, ZZ_pX b, ZZ_pX c);
inline ZZ_pX divxi(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2);
inline ZZ_pX divxi1(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2);
inline ZZ_pX divxi2(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2);
inline ZZ_pX diveta(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2);
inline ZZ_pX diveta1(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2);
inline ZZ_pX diveta2(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2);
inline ZZ_pX trunczeta(ZZ_pX a, ZZ_pX b, ZZ_pX c);
inline ZZ_pX trunczeta1(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U);
inline ZZ_pX trunczeta2(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U);
inline ZZ_pX aut(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U);
inline int autdeg(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U);
void reduce(ZZ_p half);
void reduce1(ZZ_p half, ZZ_p U);
void reduce2(ZZ_p half, ZZ_p U);
void triangle(vec_ZZ_pX A, vec_ZZ_pX &J);
int compose(vec_ZZ_pX A, vec_ZZ_pX B, vec_ZZ_pX &J);
int square(vec_ZZ_pX A, vec_ZZ_pX &J);
int relprimecomp(vec_ZZ_pX A, vec_ZZ_pX B, vec_ZZ_pX &J);
int primcomp(vec_ZZ_pX A, vec_ZZ_pX B, vec_ZZ_pX P, int g, vec_ZZ_pX &J);
int divprimes(vec_ZZ_pX &A, vec_ZZ_pX &B, ZZ_pX g, vec_ZZ_pX &I);
void babystep0(ZZ_p half, ZZ_p U);
void babystep1(ZZ_p half, ZZ_p U);
void babystep2(ZZ_p half, ZZ_p U);
int search(int lo, int hi, int verbose);
void insert(int l0, int m0, int verbose);
void reduceideal(ZZ_p half, ZZ_p U);

