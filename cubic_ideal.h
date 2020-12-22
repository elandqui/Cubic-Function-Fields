/*
  cubic_ideal.h

  -HZ
  
  Modified by Eric Landquist

*/

#ifndef CUBIC_IDEAL_H
#define CUBIC_IDEAL_H
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <iostream>

using namespace std;

NTL_CLIENT

extern ZZ_pX G;  // The defining polynomials.
extern ZZ_pX H;  // The defining polynomials.
extern ZZ_pX f;  // f = GH^2. The function field is F_q(C), C: Y^3 = f.

class cubic_ideal
{
 private:

 public:
  cubic_ideal();
  cubic_ideal(ZZ_pX& _s, ZZ_pX& _s1, ZZ_pX& _s2, ZZ_pX& _u, ZZ_pX& _v, ZZ_pX& _w);
  cubic_ideal(const cubic_ideal &A);

  //destructor
  ~cubic_ideal();

  ZZ_pX s, s1, s2, u, v, w; // basis of the ideal.
  
  //accesor functions
  // get functions
  ZZ_pX get_s() const { return s; }
  ZZ_pX get_s1() const { return s1; }
  ZZ_pX get_s2() const { return s2; } 
  ZZ_pX get_u() const { return u; }
  ZZ_pX get_v() const { return v; }
  ZZ_pX get_w() const { return w; }
  
  // set functions
  void set_s(ZZ_pX _s){ s = _s; }
  void set_s1(ZZ_pX _s1){ s1 = _s1; }
  void set_s2(ZZ_pX _s2){ s2 = _s2; } 
  void set_u(ZZ_pX _u){ u = _u; }
  void set_v(ZZ_pX _v){ v = _v; }
  void set_w(ZZ_pX _w){ w = _w; }
  
  //calculator functions
  void factors(cubic_ideal *newVal);
  
  void print();
  bool is_valid();                  // is this an honest ideal
  bool is_mult_correct(const cubic_ideal, const cubic_ideal);
  bool is_not_ramified() const ;
      
  //overloaded operators
  cubic_ideal  operator* (const cubic_ideal& I2) const;
  cubic_ideal  operator/ (const cubic_ideal& I2);
  bool operator== (const cubic_ideal &op2) const;  
  bool operator!= (const cubic_ideal &op2);
  cubic_ideal  operator^ (const ZZ &pow) const;

  void multiply (cubic_ideal &A, const cubic_ideal &B, const cubic_ideal &C);
  
  void divide (cubic_ideal &A, cubic_ideal &B, const cubic_ideal &C);
  bool IsEqual (const cubic_ideal  &A, const cubic_ideal &B);
  
  void power (cubic_ideal &A, const cubic_ideal &B, const ZZ &C);

};

// Functions for ideal reduction in unit rank 0.

void vecsort(ZZ_pX* r, long* wt);
void vecsort_general(ZZ_pX* r, long* wt);
void MinElt(ZZ_pX* min, cubic_ideal &I);
void MinElt_general(ZZ_pX* min, cubic_ideal &I);
void swap(long* a, long* b);
void swap2(long* a, long* b);
void swap(ZZ_pX* a, ZZ_pX* b);
void can_basis(cubic_ideal &I, ZZ_pX* min);
void can_basis_general(ZZ_pX &D, cubic_ideal &I, ZZ_pX* min);
void can_basis2(ZZ_pX &D, cubic_ideal &I, ZZ_pX* min);
//void can_basis(ZZ_pX &D, cubic_ideal &I, ZZ_pX* min); // Cut?
void reduce(cubic_ideal &I1, cubic_ideal &I2);

// Coefficient reduction.

void reduceCoeffs(cubic_ideal &I);

// Ideal building routines.

void ideal(cubic_ideal &I, ZZ_pX S, ZZ_pX S1, ZZ_pX U, ZZ_pX V, ZZ_pX W);
void random(cubic_ideal &A);
void naf(ZZ n, int z[], int &len);


ZZ_pX GCD3(const ZZ_pX &A, const ZZ_pX &B, const ZZ_pX &C);

// Functions for ideal multiplication.

void multV1(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void multV1a(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void multV1b(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void multV2(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void genMult(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void mult29(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void mult30(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void mult30s(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C, ZZ_pX &DD, ZZ_pX &r);
void mult31(ZZ_pX &D, ZZ_pX &DH, ZZ_pX &D1, ZZ_pX &D2, ZZ_pX &D3, ZZ_pX &D4, ZZ_pX &D5, cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void mult32(ZZ_pX &D, ZZ_pX &D1, ZZ_pX &D2, ZZ_pX &D3, cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void multiply(cubic_ideal &A, cubic_ideal  &B, cubic_ideal  &C, long &deg_diff);
void multiply_nored (cubic_ideal &A, const cubic_ideal &B, const cubic_ideal &C);

// Functions for ideal squaring.

void square(cubic_ideal &A, const cubic_ideal &B);
void square(cubic_ideal &A, const cubic_ideal &B, long &deg_diff);
void sqrV1(cubic_ideal &A, const cubic_ideal &B);
void sqrV1a(cubic_ideal &A, const cubic_ideal &B);
void sqrV1b(cubic_ideal &A, const cubic_ideal &B);
void sqrV2(cubic_ideal &A, const cubic_ideal &B);
void genSqr(cubic_ideal &A, const cubic_ideal &B);

// Functions for ideal inversion.

void inverse(cubic_ideal &A, cubic_ideal &B);
void inverse_nored(cubic_ideal &A, cubic_ideal &B);
void invV1(cubic_ideal &A, cubic_ideal &B);
void invV1a(cubic_ideal &A, cubic_ideal &B);
void invV2(cubic_ideal &A, cubic_ideal &B);
void invV2a(cubic_ideal &A, cubic_ideal &B);
void genInv(cubic_ideal &A, cubic_ideal &B);

// Functions for ideal division: only works for H=1.

void div1(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal   &C);
void divV1(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal   &C);
void divV2a(cubic_ideal  &A, const cubic_ideal   &B, const cubic_ideal  &C);
void divV2b(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void divV3(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);
void divV4(ZZ_pX &D, cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C);

// Check the order of the divisor class group of F_q(C).
int check(ZZ &order);

#endif
