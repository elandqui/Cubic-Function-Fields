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


ZZ_pX d;              // d = denominator or mu and nu
int ddeg=0;           // degree of d
ZZ_pX G;              // as above

vec_ZZ_pX I;          // C stores theta and I stores mu and nu.
vec_ZZ_pX D;          // D stores the denominators.

int p, l;             // p = preperiod length, l = period length 
int m;                // m = period lenght of second cycle

int coordl[180];      // Stores coordinates of the infrastructure.
int coordm[180];
