#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <iostream>
#include "CubicIdeal/cubic_ideal.h"
#include "Invariants/invariants.h"

using namespace std;

NTL_CLIENT

extern ZZ_pX G;   // The defining polynomials.
extern ZZ_pX H;   // The defining polynomials.
extern ZZ_pX f;   // f = GH^2. The function field is F_q(C), C: Y^3 = f.
extern int rank;  // The unit rank of F_q[C].
extern ZZ q;      // The characteristic of the field.
extern int genus; // The genus of the field.

extern int NUMPRIMES;
extern long PRIMES[];

class infrastructure_ideal
{
 private:

 public:
  infrastructure_ideal();
  infrastructure_ideal(ZZ_pX& _mu0, ZZ_pX& _mu1, ZZ_pX& _mu2, ZZ_pX& _nu0, ZZ_pX& _nu1, ZZ_pX& _nu2, ZZ_pX& d, ZZ &d0, ZZ& d1, ZZ& d2);
  infrastructure_ideal(const infrastructure_ideal &A);

  //destructor
  ~infrastructure_ideal();

  // Ideal = {1, mu0 + mu1*rho + mu2*omega, nu0 + nu1*rho + nu2*omega}/d

  ZZ_pX mu0, mu1, mu2, nu0, nu1, nu2, d; // basis of the ideal
   
  ZZ d0, d1, d2; // the distance. d1 = d2 = 0 if rank = 1.
  
  void print();
  bool operator==(infrastructure_ideal & A);
  int is_basis_reduced();
  int is_distinguished();
  // Overloaded giant step operator.
  infrastructure_ideal operator* (infrastructure_ideal& B);
};

// Initialization functions.

void init(int prec);

ZZ_pX cuberoot(ZZ_pX A);
void basis();

// Quickie functions for ideal and basis reduction. 

inline int xideg(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int xideg1(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int xideg2(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int etadeg(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int etadeg1(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int etadeg2(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int zetadeg(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int zetadeg1(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int zetadeg2(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline int discdeg(infrastructure_ideal &A);
inline int eltdeg(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline ZZ_pX divxi(infrastructure_ideal &A);
inline ZZ_pX divxi1(infrastructure_ideal &A);
inline ZZ_pX divxi2(infrastructure_ideal &A);
inline ZZ_pX diveta(infrastructure_ideal &A);
inline ZZ_pX diveta1(infrastructure_ideal &A);
inline ZZ_pX diveta2(infrastructure_ideal &A);
inline ZZ_pX trunczeta(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline ZZ_pX trunczeta1(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline ZZ_pX trunczeta2(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d);
inline ZZ_pX aut(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c);

// Computing infrastructure ideals near a given distance.

void below(ZZ &y, infrastructure_ideal &A);
void below(ZZ &y0, ZZ &y2, infrastructure_ideal &A);

// Basis reduction functions.

void reduce_basis(infrastructure_ideal &A);
void reduce_basis1(infrastructure_ideal &A);
void reduce_basis2(infrastructure_ideal &A);
void reduce_basis_nondist(infrastructure_ideal &A);
void reduce_ideal(infrastructure_ideal &A);
void reduce_ideal_r2(infrastructure_ideal &A);

// Converting an ideal from reduced to canonical bases and vice versa.

void inf_to_cubic(infrastructure_ideal &A, cubic_ideal &B);
void inf_to_cubic_gen(infrastructure_ideal &A, cubic_ideal &B);
void cubic_to_inf(cubic_ideal &A, infrastructure_ideal &B);
void cubic_to_inf_dist(cubic_ideal &A, infrastructure_ideal &B);

// Arithmetic: The Baby Step operation.

void baby_step_r1(infrastructure_ideal &A, infrastructure_ideal &B);
void baby_step_0_r2(infrastructure_ideal &A, infrastructure_ideal &B);
void baby_step_1_r2(infrastructure_ideal &A, infrastructure_ideal &B);
void baby_step_2_r2(infrastructure_ideal &A, infrastructure_ideal &B);

// Arithmetic: The Giant Step operation.

void giant_step_r1(infrastructure_ideal &result, infrastructure_ideal &A, cubic_ideal& B, ZZ &Bdist);
void giant_step_r2(infrastructure_ideal &result, infrastructure_ideal &A, cubic_ideal& B, ZZ &Bd0, ZZ &Bd1, ZZ &Bd2);

// Arithmetic: The Inverse operation.

void inverse(infrastructure_ideal &A, infrastructure_ideal &B);

// Functions for extracting the regulator from a multiple.
ZZ extract(ZZ h0, ZZ bound, vec_ZZ& factors, vec_long& exponents, int manext);
ZZ extract2(ZZ h0, ZZ bound, vec_ZZ& factors, vec_long& exponents, int manext);
int factor(ZZ n, vec_ZZ& factors, vec_long& exponents);
void pollard_rho(ZZ n, ZZ& prime);

// Check if h is a multiple of the regulator.
int check_inf(ZZ &h);
int check_inf2(ZZ &h);
