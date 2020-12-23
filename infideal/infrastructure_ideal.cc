/* 
   infrastructure_ideal.cc

   This file contains functions to perform infrastructure arithmetic such as 
   baby steps, giant steps, and inversion. Also included are routines to 
   reduce ideals and bases, along with basis conversion functions and a
   method to extract the regulator from a multiple. 
   For the functions, we will use results and algorithms from:

   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)
   [LSY] Y. Lee, R. Scheidler, and C. Yarrish, "Computation of the 
                        fundamental units and the regulator of a cyclic 
			cubic function field"
   [S01] R. Scheidler, "Ideal arithmetic and infrastructure in purely 
                        cubic function fields"
   
   How to create the library:
   --------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download cubic_ideal.cpp, cubic_ideal.h and makefile into 
      ../CubicIdeal.
   Compile libcubic via make in ../CubicIdeal.
   Download infrastructure_ideal.cc, infrastructure_ideal.h and makefile 
      into ../Infideal.
   Compile libinfideal via make in ../Infideal.

   Make sure that the paths to the NTL libraries are
   specified correctly in the makefile.
   
   Known Problems: 
   ---------------
   In below() for unit rank 2, there is an upper bound on y2 beyond which
   the function produces incorrect output (or perhaps none at all).
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   4. Fix below() to work more generally in unit rank 2.

   Authors: 
   Eric Landquist
   Christopher Yarrish 
*/

#include "infrastructure_ideal.h"

NTL_CLIENT

ZZ_pX rho;       // rho = (x^2G(x))^(1/3): basis element of F_q[C]
ZZ_pX omega;     // omega = (xG^2(x))^(1/3): basis element of F_q[C]
int prec;        // The precision of the Puiseux series approximations.
int belowinit;
ZZ UU;           // A primitive cube root of 1 mod q.
ZZ UU2;          // UU^2 mod q (i.e. the other primitive cube root of 1.)
ZZ half;         // 2^{-1} mod q  
ZZ One = to_ZZ(1);
ZZ Zero = ZZ();

//default constructor
infrastructure_ideal::infrastructure_ideal() {
  mu0 = 0; mu1 = 1;  mu2 = 0;  nu0 = 0;  nu1 = 0;  nu2 = 1;  d = 1;
  d0 = d1 = d2 = 0;
}

// constructor
infrastructure_ideal::infrastructure_ideal(ZZ_pX& _mu0, ZZ_pX& _mu1, ZZ_pX& _mu2, ZZ_pX& _nu0, ZZ_pX& _nu1, ZZ_pX& _nu2, ZZ_pX& _d, ZZ& _d0, ZZ& _d1, ZZ& _d2){
  mu0 = _mu0;  mu1 = _mu1;  mu2 = _mu2;  nu0 = _nu0;  nu1 = _nu1;  nu2 = _nu2;  d = _d;  
  d0 = _d0;  
  d1 = _d1;
  d2 = _d2;
}

// constructor
infrastructure_ideal::infrastructure_ideal(const infrastructure_ideal &A){
  mu0 = A.mu0;  mu1 = A.mu1;  mu2 = A.mu2;  nu0 = A.nu0;  nu1 = A.nu1;  nu2 = A.nu2;  d = A.d;
  d0 = A.d0;
  d1 = A.d1;
  d2 = A.d2;
}

//destructor
infrastructure_ideal::~infrastructure_ideal() {

}

bool infrastructure_ideal::operator==(infrastructure_ideal  & A){
  if((this->mu0 == A.mu0) && (this->mu1 == A.mu1) && (this->mu2 == A.mu2) && (this->nu0 == A.nu0) && (this->nu1 == A.nu1) && (this->nu2 == A.nu2) && (this->d == A.d))
    return true;
  else
    return false;
}

void infrastructure_ideal::print(){
  cout<<this->mu0<<endl;
  cout<<this->mu1<<endl;
  cout<<this->mu2<<endl;
  cout<<this->nu0<<endl;
  cout<<this->nu1<<endl;
  cout<<this->nu2<<endl;
  cout<<this->d<<endl;
}

// Input: an infrastructure ideal.
// Output: 0 if the basis of the ideal is not 0-reduced.
//         1 if the basis of the ideal is 0-reduced.
// See Eqn. 4-2 of [LSY] or Def. 5.1.1 of [L09].

int infrastructure_ideal::is_basis_reduced(){
  
  //cout<<xideg(this->nu1, this->nu2, this->d)<<" < "<<xideg(this->mu1, this->mu2, this->d)<<endl;
  //cout<<etadeg(this->mu1, this->mu2, this->d)<<" < 0 <= "<<etadeg(this->nu1, this->nu2, this->d)<<endl;
  //cout<<zetadeg(this->nu0, this->nu1, this->nu2, this->d)<<" < 0 > "<<zetadeg(this->mu0, this->mu1, this->mu2, this->d)<<endl;
  //cout<<"Is ideal distinguished? "<<eltdeg(this->mu0, this->mu1, this->mu2, this->d)<<" > 0 ? max "<<eltdeg(this->nu0, this->nu1, this->nu2, this->d)<<", "<<etadeg(this->nu1, this->nu2, this->d)<<" > 0?"<<endl;
  
  if(xideg(this->nu1, this->nu2, this->d) >= xideg(this->mu1, this->mu2, this->d)){
    //cout<<"xideg: "<<xideg(this->nu1, this->nu2, this->d)<<" < "<<xideg(this->mu1, this->mu2, this->d)<<endl;
    return 0;
  }
  if(etadeg(this->mu1, this->mu2, this->d) >= 0){
    //cout<<"etamudeg: "<<etadeg(this->mu1, this->mu2, this->d)<<" < 0"<<endl;
    return 0;
  }
  if(etadeg(this->nu1, this->nu2, this->d) < 0){ 
    //cout<<"etanudeg: 0 <= "<<etadeg(this->nu1, this->nu2, this->d)<<endl;
    return 0;
  }
  if(zetadeg(this->mu0, this->mu1, this->mu2, this->d) >= 0){ 
    //cout<<"zetamudeg: "<<zetadeg(this->mu0, this->mu1, this->mu2, this->d)<<" < 0"<<endl;
    return 0;
  }
  if(zetadeg(this->nu0, this->nu1, this->nu2, this->d) > 0){
    //cout<<"zetanudeg: "<<zetadeg(this->nu0, this->nu1, this->nu2, this->d)<<" < 0"<<endl;
    return 0;
  }
  if(etadeg(this->nu1, this->nu2, this->d) == 0){
    if(eltdeg(this->nu0, this->nu1, this->nu2, this->d) == 0){
      //cout<<"max "<<eltdeg(this->nu0, this->nu1, this->nu2, this->d)<<", "<<etadeg(this->nu1, this->nu2, this->d)<<" > 0?"<<endl;
      return 0;
    }
  }
  return 1;
}

// Input: an infrastructure ideal with a 0-reduced basis.
// Output: 0 if the ideal is not distinguished.
//         1 if the ideal is distinguished. (i.e. an infrastructure ideal)
// See Thm 5.2.1 of [L09] (Thm. 6.6 of [S01] for rank = 1)

int infrastructure_ideal::is_distinguished(){
  if(is_basis_reduced()){
    if(eltdeg(this->mu0, this->mu1, this->mu2, this->d) <= 0){
      //cout<<"Is ideal distinguished? "<<eltdeg(this->mu0, this->mu1, this->mu2, this->d)<<" > 0 ? max "<<eltdeg(this->nu0, this->nu1, this->nu2, this->d)<<", "<<etadeg(this->nu1, this->nu2, this->d)<<" > 0?"<<endl;
      return 0;
    }
    else if( (eltdeg(this->nu0, this->nu1, this->nu2, this->d) > 0) || (etadeg(this->nu1, this->nu2, this->d) > 0) ){
      return 1;
    }
    else{
      //cout<<"Is ideal distinguished? "<<eltdeg(this->mu0, this->mu1, this->mu2, this->d)<<" > 0 ? max "<<eltdeg(this->nu0, this->nu1, this->nu2, this->d)<<", "<<etadeg(this->nu1, this->nu2, this->d)<<" > 0?"<<endl;
      return 0;
    }
  }
  else{
    //cout<<"Is ideal distinguished? "<<eltdeg(this->mu0, this->mu1, this->mu2, this->d)<<" > 0 ? max "<<eltdeg(this->nu0, this->nu1, this->nu2, this->d)<<", "<<etadeg(this->nu1, this->nu2, this->d)<<" > 0?"<<endl;
    return 0;
  }
  return 0;
}

/***********************************************/
/* Initialize the infrastructure with:         */
/*   1. precision                              */
/*   2. computing the basis of O = Fq[C]       */
/*   3. optimizing the running time of below() */
/*   4. If rank=2, find a cube root of 1       */
/*                 and compute 2^{-1} mod q.   */
/***********************************************/
void init(int _prec){
  prec = _prec;
  basis();
  belowinit = belowsteps();
  ZZ_p _UU, _half;

  if(rank == 2){
    _UU = to_ZZ_p(1);
    ZZ_p b = to_ZZ_p(1);
    ZZ third = (q-1)/3;
    while(IsOne(_UU)){
      b++;
      power(_UU, b, third);
    }
    UU = rep(_UU);
    UU2 = sqr(UU)%q;

    if(UU2 < UU){
      UU += UU2;
      UU2 = UU - UU2;
      UU -= UU2;
    }
  }
    _half = inv(to_ZZ_p(2));
    half = rep(_half);
}

// Compute an approximation of A^(1/3) with prec precision.
// For computing the basis of F_q[C].
// Written by Christopher Yarrish in Maple, converted to C++ and 
//  modified by EJL.

ZZ_pX cuberoot(ZZ_pX A){
  ZZ_pX C;         // the cube root
  ZZ_pX T = A;     // temp polynomial
  int i, j, r, s, t;
  ZZ_p temp, denom;
  int k = deg(A);
  vec_ZZ_p Dseq; 
  vec_ZZ_p cube;
  int cubedeg = k/3;
  int ind = cubedeg + prec+1;

  cube.SetLength(ind);
  cube[ind-1] = 1;

  Dseq.SetLength(k+1);

  for(i=0; i<=k; i++){
    Dseq[k-i] = coeff(T, k-i);
    if( !IsZero(Dseq[k-i]) )
      T-=ZZ_pX(deg(T), LeadCoeff(T));
    if( IsZero(coeff(A, 0)) )
      Dseq[0] = 0;
  }

  for(i=1; i<ind; i++){
    // First term
    if(k-i >= 0)
      temp = Dseq[k-i];
    else
      temp=0;

    // Second term
    for(t=ind-i+3; t<=ind; t++){
      for(s=ind-i+2; s<t; s++){
	for(r=ind-i+1; r<s; r++){
	  if(((r+s+t)-(3*prec+3)) == (k-i))
	    temp -= 6*cube[r-1]*cube[s-1]*cube[t-1];
	}
      }
    }

    // Third term
    for(r=ind-i+1; r<=ind; r++){
      for(s=ind-i+1; s<=ind; s++){
	if(((2*r+s)-3*(prec+1) == (k-i)) && (r!=s))
	  temp -= 3*sqr(cube[r-1])*cube[s-1];
      }
    }

    // Fourth term
    if(!(i%3))
      temp -= power(cube[cubedeg+prec-(i/3)], 3);

    //Finally
    if(i==1){
      denom = 3*sqr(cube[ind-1]);
      denom = inv(denom);
    }
    temp*=denom;

    cube[ind-i-1] = temp;
  }
  for(j=0; j<ind; j++)
    C+=ZZ_pX(j, cube[j]);

  return C;
}

// Computes the basis, [1, rho, omega], of F_q[C].

void basis(){

  rho = cuberoot(f);
  omega = cuberoot(sqr(G)*H);

}

/*

  The following several functions operate on function field elements,
  alpha = a + b*rho + c*omega \in F_q(C) needed for basis and ideal
  reduction.
  See Eqn. 4-1 of [LSY] or Sect. 5.1.1 of [L09] for the definitions
  of xi, eta, and zeta.

 */

// xi_alpha = b*rho + c*omega
inline int xideg(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  return(deg(b*rho + c*omega) - prec - deg(d));
}

inline int xideg1(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  return(deg(b*to_ZZ_p(UU)*rho + c*to_ZZ_p(UU2)*omega) - prec - deg(d));
}

inline int xideg2(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
    return(deg(b*to_ZZ_p(UU2)*rho + c*to_ZZ_p(UU)*omega) - prec - deg(d));
}

// eta_alpha = b*rho - c*omega
inline int etadeg(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
//cout<<"eta: "<<b*rho - c*omega<<endl;
  return(deg(b*rho - c*omega) - prec - deg(d));
}

inline int etadeg1(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  //cout<<"eta: "<<b*rho - c*omega<<endl;
  return(deg(b*to_ZZ_p(UU)*rho - c*to_ZZ_p(UU2)*omega) - prec - deg(d));
}

inline int etadeg2(ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  //cout<<"eta: "<<b*rho - c*omega<<endl;
  return(deg(b*to_ZZ_p(UU2)*rho - c*to_ZZ_p(UU)*omega) - prec - deg(d));
}


inline int zetadeg(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  return(deg(LeftShift(2*a, prec) - b*rho - c*omega) - prec - deg(d));
}

inline int zetadeg1(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  return(deg(LeftShift(2*a, prec) - b*to_ZZ_p(UU)*rho - c*to_ZZ_p(UU2)*omega) - prec - deg(d));
}

inline int zetadeg2(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  return(deg(LeftShift(2*a, prec) - b*to_ZZ_p(UU2)*rho - c*to_ZZ_p(UU)*omega) - prec - deg(d));
}

// The degree of the discriminant of the ideal
inline int discdeg(infrastructure_ideal &A){
  return (2*(deg((A.mu1*rho + A.mu2*omega)*(A.nu1*rho - A.nu2*omega) - (A.mu1*rho - A.mu2*omega)*(A.nu1*rho + A.nu2*omega)) - 2*prec));
}

inline int eltdeg(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  return(deg(LeftShift(a, prec) + b*rho + c*omega) - prec - deg(d));
}

inline ZZ_pX divxi(infrastructure_ideal &A){
  ZZ_pX ximu, xinu, k;
  ximu = A.mu1*rho + A.mu2*omega;
  xinu = A.nu1*rho + A.nu2*omega;

  div(k, ximu, xinu);
  return k;
}

inline ZZ_pX divxi1(infrastructure_ideal &A){
  ZZ_pX ximu, xinu, k;
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);

  ximu = A.mu1*_UU*rho + A.mu2*_UU2*omega;
  xinu = A.nu1*_UU*rho + A.nu2*_UU2*omega;

  div(k, ximu, xinu);
  return k;
}

inline ZZ_pX divxi2(infrastructure_ideal &A){
  ZZ_pX ximu, xinu, k;
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);

  ximu = A.mu1*_UU2*rho + A.mu2*_UU*omega;
  xinu = A.nu1*_UU2*rho + A.nu2*_UU*omega;

  div(k, ximu, xinu);
  return k;
}

inline ZZ_pX diveta(infrastructure_ideal &A){
  ZZ_pX etamu, etanu, k;
  
  etamu = A.mu1*rho - A.mu2*omega;
  etanu = A.nu1*rho - A.nu2*omega;
  
  div(k, etanu, etamu);

  return k;
}

inline ZZ_pX diveta1(infrastructure_ideal &A){
  ZZ_pX etamu, etanu, k;
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);

  etamu = A.mu1*_UU*rho - A.mu2*_UU2*omega;
  etanu = A.nu1*_UU*rho - A.nu2*_UU2*omega;
 
  div(k, etanu, etamu);

  return k;
}

inline ZZ_pX diveta2(infrastructure_ideal &A){
  ZZ_pX etamu, etanu, k;
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);

  etamu = A.mu1*_UU2*rho - A.mu2*_UU*omega;
  etanu = A.nu1*_UU2*rho - A.nu2*_UU*omega;
 
  div(k, etanu, etamu);

  return k;
}

inline ZZ_pX trunczeta(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  ZZ_pX k;
  
  div(k, LeftShift(2*a, prec) - b*rho - c*omega, LeftShift(d, prec));

  return(k*d);
}

inline ZZ_pX trunczeta1(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  ZZ_pX k;
  
  div(k, LeftShift(2*a, prec) - b*to_ZZ_p(UU)*rho - c*to_ZZ_p(UU2)*omega, LeftShift(d, prec));

  return(k*d);
}

inline ZZ_pX trunczeta2(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c, ZZ_pX &d){
  ZZ_pX k;
  
  div(k, LeftShift(2*a, prec) - b*to_ZZ_p(UU2)*rho - c*to_ZZ_p(UU)*omega, LeftShift(d, prec));

  return(k*d);
}

inline ZZ_pX aut(ZZ_pX &a, ZZ_pX &b, ZZ_pX &c){
  return(a + RightShift(to_ZZ_p(UU)*b*rho + to_ZZ_p(UU2)*c*omega, prec));
}

// Algorithm 5.3.26 of [Landquist, 2009]
// Input: an integer y, unit rank 1.
// Output: the ideal E = D(y), where E.d0 <= y < babystep(E).d0.

void below(ZZ &y, infrastructure_ideal &E){
  ZZ s, t, cap;
  infrastructure_ideal A, D, D0, Inv;
  int i, len;
  
  E.mu0 = 0; E.mu1 = 1;  E.mu2 = 0;  
  E.nu0 = 0;  E.nu1 = 0;  E.nu2 = 1;  
  E.d = 1;
  E.d0 = 0;
  
  //E.print();
  reduce_basis(E);
  //E.print();
  for(i=0; i< belowinit; i++){
    baby_step_r1(E, D0);
    E = D0;
  }
  //D0.print();
  s = y/D0.d0;
  //cout<<y<<" "<<D0.d0<<" "<<s<<" "<<Inv.d0<<endl;
  if(IsZero(s)){
    E.mu0 = 0; E.mu1 = 1;  E.mu2 = 0;  
    E.nu0 = 0;  E.nu1 = 0;  E.nu2 = 1;  
    E.d = 1;
    E.d0 = 0;
    reduce_basis(E);
    if(E.d0 <= y){
      while(E.d0 <= y){
	D0 = E;
	baby_step_r1(D0, E);
      }
      E = D0;
    }
    return;
  }

  inverse(Inv, D0);

  //cout<<y<<" "<<D0.d0<<" "<<s<<" "<<Inv.d0<<endl;

  // To store NAF(s). 
  int z[NumBits(s)+2];
  naf(s, z, len);
  //for(int j=len-1; j>=0; j--)
  //  cout<<z[j]<<" ";
  //cout<<endl;
  t=z[len-1];
  
  for(i=len-2; i>=0; i--){
    A = E;
    E = A*A;
    cap = 2*t*D0.d0;
    //cout<<E.d0<<" "<<E.d<<" spot 1"<<endl;
    if(E.d0 <= cap){
      while(E.d0 <= cap){
	D = E;
	baby_step_r1(D, E);
      }
      E = D;
    }
    t = 2*t + z[i];
    if(z[i]==1){
      A = E;
      E = A*D0;
      //cout<<"A:  "<<A.d0<<" "<<A.d<<endl;
      //cout<<"D0: "<<D0.d0<<" "<<D0.d<<endl;
      //cout<<E.d0<<" "<<E.d<<" spot 2"<<endl;
      cap = t*D0.d0; 
      if(E.d0 <= cap){
	while(E.d0 <= cap){
	  D = E;
	  baby_step_r1(D, E);
	}
	E=D;
      }
    }
    else if(z[i]==-1){
      A = E;
      E = A*Inv;
      //cout<<E.d0<<" "<<E.d<<" spot 3"<<endl;
      cap = t*D0.d0;
      //cout<<E.d0<<" "<<D0.d0<<" "<<t<<" "<<cap<<" "<<A.d0<<" inverse"<<endl;
      while(E.d0 > cap){
	A = E;
	E = A*Inv;
	//cout<<E.d0<<" "<<E.d<<" spot 4"<<endl;
      }
      //cout<<E.d0<<" "<<D0.d0<<" "<<t<<" "<<cap<<" "<<A.d0<<" inverse"<<endl;
      if(E.d0 <= cap){
	while(E.d0 <= cap){
	  D = E;
	  baby_step_r1(D, E);
	}
	E=D;
      }
    }
    else;
  }
  
  while(E.d0 <= y){
    if(E.d0 == y){
      return;
    }
    D = E;
    baby_step_r1(D, E);
  }

  E = D;
  return;
}

// Generalization of Algorithm 5.3.26 of [Landquist, 2009] to r = 2.
// Input: integers y0, y2.
// Output: an ideal E = D(y0, y2), where 
//         E.d0 <= y0 < bs_0(E).d0 and E.d2 <= y2 < bs_2(E).d2.

void below(ZZ &y0, ZZ &y2, infrastructure_ideal &E){
  ZZ s, t, cap;
  infrastructure_ideal A, B, D, D0, Inv;
  int i, len;//, j;

  E.mu0 = 0; E.mu1 = 1;  E.mu2 = 0;  
  E.nu0 = 0;  E.nu1 = 0;  E.nu2 = 1;  
  E.d = 1;
  E.d0 = E.d1 = E.d2 = 0;

  // First find below(y0, 0).

  reduce_basis(E);

  for(i=0; i< belowinit; i++){
    baby_step_0_r2(E, D0);
    E = D0;
  }
  
  if(E.d2 < 0){
    reduce_basis2(E);
    while(E.d2 <= 0){   
      B=E;
      baby_step_2_r2(B, E);
    }
    reduce_basis(B);
    E = B;
    if(E.d2 < 0){
      do{
	baby_step_0_r2(E, B);
	E = B;
	reduce_basis2(B);
	while(B.d2 <= 0){
	  A = B;
	  baby_step_2_r2(A, B);
	}
	E = B;
	reduce_basis(E);
      } while(E.d2 < 0);
    }
  }
  //cout<<"A: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  s = y0/E.d0;
  
  if(IsZero(s)){
    E.mu0 = 0; E.mu1 = 1;  E.mu2 = 0;  
    E.nu0 = 0;  E.nu1 = 0;  E.nu2 = 1;  
    E.d = 1;
    E.d0 = 0;
    reduce_basis(E);
    if(E.d0 <= y0){
      while(E.d0 <= y0){
	D0 = E;
	baby_step_0_r2(D0, E);
      }
      E = D0;
    }
    if(E.d2 < 0){
      reduce_basis2(E);
      while(E.d2 <= 0){   
	B=E;
	baby_step_2_r2(B, E);
    }
      reduce_basis(B);
      E = B;
      if(E.d2 < 0){
	do{
	  baby_step_0_r2(E, B);
	  E = B;
	  reduce_basis2(B);
	  while(B.d2 <= 0){
	    A = B;
	    baby_step_2_r2(A, B);
	  }
	  E = B;
	  reduce_basis(E);
	} while(E.d2 < 0);
      }
    }
    
    return;
  }
  
  inverse(Inv, D0);

  // To store NAF(s). 
  int z[NumBits(s)+2];
  naf(s, z, len);
  t=z[len-1];
  
  for(i=len-2; i>=0; i--){
    A = E;
    E = A*A;
 
    cap = 2*t*D0.d0;
    
    if(E.d0 <= cap){
      while(E.d0 <= cap){
	D = E;
	baby_step_0_r2(D, E);
      }
      E = D;
    }
    
    if(E.d2 < 0){
      reduce_basis2(E);
      while(E.d2 <= 0){   
	B=E;
	baby_step_2_r2(B, E);
      }
      reduce_basis(B);
      E = B;
      if(E.d2 < 0){
	do{
	  baby_step_0_r2(E, B);
	  E = B;
	  reduce_basis2(B);
	  while(B.d2 <= 0){
	  A = B;
	  baby_step_2_r2(A, B);
	  }
	  E = B;
	  reduce_basis(E);
	} while(E.d2 < 0);
      }
    }
    //cout<<"B: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    t = 2*t + z[i];
    if(z[i]==1){
      A = E;
      E = A*D0;
      
      cap = t*D0.d0;

      if(E.d0 <= cap){
	while(E.d0 <= cap){
	  D = E;
	  baby_step_0_r2(D, E);
	}
	E = D;
      }

      if(E.d2 < 0){
	reduce_basis2(E);
	while(E.d2 <= 0){   
	  B=E;
	  baby_step_2_r2(B, E);
	}
	reduce_basis(B);
	E = B;
	if(E.d2 < 0){
	  do{
	    baby_step_0_r2(E, B);
	    E = B;
	    reduce_basis2(B);
	    while(B.d2 <= 0){
	      A = B;
	      baby_step_2_r2(A, B);
	    }
	    E = B;
	    reduce_basis(E);
	  } while(E.d2 < 0);
	}
      }   
      //cout<<"C: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;  
    }
    else if(z[i]==-1){
      A = E;
      E = A*Inv;
      
      cap = t*D0.d0;
     
      while(E.d0 > cap){
	A = E;
	E = A*Inv;
      }
      
      while(E.d0 <= cap){
	D = E;
	baby_step_0_r2(D, E);
      }
      E=D;
      //cout<<"d: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
      if(E.d2 > 0){
	reduce_basis1(E);
	while(E.d2 >= 0){
	  B=E;
	  baby_step_1_r2(B, E);
	  //cout<<"dd: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
	}
	if(B.d2 == 0)
	  E = B;
	reduce_basis(E);
      }
      else if(E.d2 < 0){
	reduce_basis2(E);
	while(E.d2 <= 0){   
	  B=E;
	  baby_step_2_r2(B, E);
	}
	reduce_basis(B);
	E = B;
	if(E.d2 < 0){
	  do{
	    baby_step_0_r2(E, B);
	    E = B;
	    reduce_basis2(B);
	    while(B.d2 <= 0){
	      A = B;
	      baby_step_2_r2(A, B);
	    }
	    E = B;
	    reduce_basis(E);
	  } while(E.d2 < 0);
	}
      }     
      else;
      //cout<<"D: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    }
    else;
  }
  //cout<<"1: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  if(E.d0 < y0){
    while(E.d0 <= y0){
      if(E.d0 == y0 && E.d2 == y2){
	return;
      }
      D = E;
      baby_step_0_r2(D, E);
      //cout<<"1a: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    }
    if((D.d2 < 0) && (E.d2 < 0)){
      reduce_basis2(E);
      while(E.d2 <= 0){
	D = E;
	baby_step_2_r2(D, E);
      }
      reduce_basis(D);
      E = D;
    }
    //cout<<"1b: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    E = D;
  }
  //cout<<"2: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  if(E.d2 < 0){
    reduce_basis2(E);
    while(E.d2 <= 0){   
      B=E;
      baby_step_2_r2(B, E);
    }
    reduce_basis(B);
    E = B;
    if(E.d2 < 0){
      do{
	baby_step_0_r2(E, B);
	E = B;
	reduce_basis2(B);
	while(B.d2 <= 0){
	  A = B;
	  baby_step_2_r2(A, B);
	}
	E = B;
	reduce_basis(E);
      } while(E.d2 < 0);
    }
  }    
  //cout<<"3: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl; 
  
  if(IsZero(y2))
    return;
  //cout<<"After first stage of below."<<endl;
  //cout<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  //E.print();

  // Now we proceed to find the y2 part.
  
  if(!IsZero(y2)){// <= 20000){
    reduce_basis2(E);
    while(E.d2 <= y2){
      D0 = E;
      //reduce_basis2(E);
      baby_step_2_r2(D0, E);
      //cout<<D0.d0<<" "<<D0.d1<<" "<<D0.d2<<endl;
      //cout<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
      //reduce_basis(D0);
      //baby_step_0_r2(D0, E);
    }
    E = D0;
    //cout<<"Got the d2"<<endl;
    //cout<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    if(E.d0 < y0){
      reduce_basis(E);
      while(E.d0 <= y0){
	D0 = E;
	baby_step_0_r2(D0, E);
      }
      if((D0.d2 < y2) && (E.d2 < y2)){
	reduce_basis2(E);
	while(E.d2 <= y2){
	  D0 = E;
	  baby_step_2_r2(D0, E);
	}
	reduce_basis(D0);
	E = D0;
      }
      E = D0;
    }
    //cout<<"1: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<" "<<(E.d0 == y0)<<endl;
    for(int k=0; (k<5) && !((E.d0 == y0) && (E.d2 == y2)) ; k++){
      if(E.d2 < y2){
	reduce_basis2(E);
	while(E.d2 <= y2){
	  D0 = E;
	  baby_step_2_r2(D0, E);
	  //cout<<"2: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
	}
	//cout<<"3: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
	E = D0;
	//cout<<"4: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
	//baby_step_0_r2(D0, E);
      }
      if(E.d0 < y0){
	reduce_basis(E);
	while(E.d0 <= y0){
	  D0 = E;
	  baby_step_0_r2(D0, E);
	  //cout<<"5: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
	}
	if((D0.d2 < y2) && (E.d2 < y2)){
	  reduce_basis2(E);
	  while(E.d2 <= y2){
	    D0 = E;
	    baby_step_2_r2(D0, E);
	  }
	  reduce_basis(D0);
	  E = D0;
	}
	E = D0;
	//baby_step_0_r2(D0, E);
      }
    }
    return; 
  }

  // BELOW NEEDS FIXING
  /*
  //E.print();
  for(i=0; i< belowinit; i++){
    reduce_basis2(E);
    baby_step_2_r2(E, D0);
    reduce_basis(D0);
    baby_step_0_r2(D0, E);
  }
  //cout<<"1: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  if(E.d0 < y0){
    reduce_basis(E);
    while(E.d0 <= y0){
      D0=E;
      baby_step_0_r2(D0, E);
    }
  }
  D0 = E;
  //cout<<"2: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  //D0.print();
  s = y2/D0.d2;
  //cout<<y<<" "<<D0.d0<<" "<<s<<" "<<Inv.d0<<endl;
  //cout<<"3: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  inverse(Inv, D0);

  //cout<<y<<" "<<D0.d0<<" "<<s<<" "<<Inv.d0<<endl;

  // To store NAF(s). 
  int z2[NumBits(s)+2];
  naf(s, z2, len);
  //for(int j=len-1; j>=0; j--)
  //  cout<<z2[j]<<" ";
  //cout<<endl;
  t=z2[len-1];
  
  for(i=len-2; i>=0; i--){
    A = E;
    E = A*A;
 
    cap = 2*t*D0.d2;
    //cout<<E.d0<<" "<<E.d<<" spot 1"<<endl;
    while(E.d2 <= cap){
      reduce_basis2(E);
      baby_step_2_r2(E, D);
      reduce_basis(D);
      baby_step_0_r2(D, E);
    }
    //cout<<"4: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    if(E.d0 < y0){
      reduce_basis(E);
      while(E.d0 <= y0){
	B=E;
	reduce_basis(B);
	baby_step_0_r2(B, E);
      }
      E = B;
    }
    //cout<<"5: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    t = 2*t + z2[i];
    if(z2[i]==1){
      A = E;
      E = A*D0;
      //cout<<"A:  "<<A.d0<<" "<<A.d<<endl;
      //cout<<"D0: "<<D0.d0<<" "<<D0.d<<endl;
      //cout<<E.d0<<" "<<E.d<<" spot 2"<<endl;
      cap = t*D0.d2; 
      while(E.d2 <= cap){
	reduce_basis2(E);
	baby_step_2_r2(E, D);
	reduce_basis(D);
	baby_step_0_r2(D, E);
      }
      //cout<<"6: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
      if(E.d2 < 0){
	reduce_basis2(E);
	while(E.d2 <= 0){
	  B=E;
	  baby_step_2_r2(B, E);
	  count++;
	}
	reduce_basis(B);
	for(j=0; j<count-1; j++){
	  E = B;
	  baby_step_0_r2(E, B);
	}
	E = B;
	count=0;
      }
      //cout<<"6: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    }
    else if(z2[i]==-1){
      A = E;
      E = A*Inv;
      //cout<<E.d0<<" "<<E.d<<" spot 3"<<endl;
      cap = t*D0.d0;
      //cout<<E.d0<<" "<<D0.d0<<" "<<t<<" "<<cap<<" "<<A.d0<<" inverse"<<endl;
      while(E.d0 > cap){
	A = E;
	E = A*Inv;
	//cout<<E.d0<<" "<<E.d<<" spot 4"<<endl;
      }
      //cout<<E.d0<<" "<<D0.d0<<" "<<t<<" "<<cap<<" "<<A.d0<<" inverse"<<endl;
      while(E.d0 <= cap){
	D = E;
	baby_step_0_r2(D, E);
      }
      E=D;
      //cout<<"7: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
      if(E.d2 > 0){
	reduce_basis1(E);
	while(E.d2 >= 0){
	  B=E;
	  baby_step_1_r2(B, E);
	}
	reduce_basis(B);
	E = B;
      }
      else if(E.d2 < 0){
	reduce_basis2(E);
	while(E.d2 <= 0){
	  B=E;
	  baby_step_2_r2(B, E);
	  count++;
	}
	reduce_basis(B);
	for(j=0; j<count-1; j++){
	  E = B;
	  baby_step_0_r2(E, B);
	}
	E = B;
	count=0;
      }
      else;
      //cout<<"8: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
    }
    else;
  }
  
  while(E.d0 <= y0){
    if(E.d0 == y0 && E.d2 == y2){
      return;
    }
    D = E;
    baby_step_0_r2(D, E);
  }
  E = D;
  //cout<<"9: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  if(E.d2 < 0){
    reduce_basis2(E);
    while(E.d2 <= 0){
      B=E;
      baby_step_2_r2(B, E);
      count++;
    }
    reduce_basis(B);
    for(j=0; j<count-1; j++){
      E = B;
      baby_step_0_r2(E, B);
    }
    E = B;
    count=0;
  }
  */
  //cout<<"10: "<<E.d0<<" "<<E.d1<<" "<<E.d2<<endl;
  return;
}

// Input: A distinguished infrastructure ideal, A. 
// Output: The ideal A given with the (0-)reduced basis.
// See Alg. 4.6 of [LSY] or Alg. 5.1.5 of [L09].

void reduce_basis(infrastructure_ideal &A){
  int xm, xn, em, en; 
  ZZ_pX temp0, temp1, temp2, k, etamu, etanu;
  ZZ_p lc, lc1;
  ZZ_p _half = to_ZZ_p(half);

  //Step 2 in 4.6
  xm = xideg(A.mu1, A.mu2, A.d);
  xn = xideg(A.nu1, A.nu2, A.d);

  if(xm < xn){
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0;  A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0 = -temp0; A.nu1 = -temp1; A.nu2 = -temp2;
  }
  if(xm == xn){
    em = etadeg(A.mu1, A.mu2, A.d); en = etadeg(A.nu1, A.nu2, A.d);
    if(em < en){
      temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
      A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
      A.nu0 = -temp0; A.nu1 = -temp1; A.nu2 = -temp2;
    }
  }

  // Step 3.
  em = etadeg(A.mu1, A.mu2, A.d); en = etadeg(A.nu1, A.nu2, A.d);
  if(em >= en){
    // Step 3.1
    while( deg(sqr(A.nu1*rho) - sqr(A.nu2*omega)) - 2*prec  > (double)discdeg(A)/2  ){
      k = divxi(A);
      temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
      A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
      A.nu0 = k*A.mu0-temp0; A.nu1 = k*A.mu1-temp1; A.nu2 = k*A.mu2-temp2;  
    }
    // Step 3.2
    k = divxi(A);
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0 = k*A.mu0-temp0; A.nu1 = k*A.mu1-temp1; A.nu2 = k*A.mu2-temp2;

    // Step 3.3
    etamu = A.mu1*rho - A.mu2*omega;
    etanu = A.nu1*rho - A.nu2*omega;
    
    em = deg(etamu); en = deg(etanu);
    if(em==en){
      lc=ZZ_p(LeadCoeff(etamu));
      lc1=ZZ_p(LeadCoeff(etanu));
      lc*=inv(lc1);
      
      A.mu0 -= lc1*A.nu0; A.mu1 -= lc1*A.nu1; A.mu2 -= lc1*A.nu2; 
    }
  }
  //Step 4.
  
  while(etadeg(A.mu1, A.mu2, A.d) >= 0){
    k = diveta(A);
    temp0 = A.nu0; temp1 = A.nu1; temp2 = A.nu2;
    A.nu0 = A.mu0; A.nu1 = A.mu1; A.nu2 = A.mu2;
    A.mu0*=k; A.mu1*=k; A.mu2*=k;
    A.mu2-=temp2; A.mu1-=temp1; A.mu0-=temp0;
  }
  // Step 5.
  if(zetadeg(A.mu0, A.mu1, A.mu2, A.d) >= 0)
    A.mu0 -= _half*trunczeta(A.mu0, A.mu1, A.mu2, A.d);
  if(zetadeg(A.nu0, A.nu1, A.nu2, A.d) >= 0)
    A.nu0 -= _half*trunczeta(A.nu0, A.nu1, A.nu2, A.d);
  
  // Normalize so mu and nu are monic.
  lc = ZZ_p(LeadCoeff(A.mu0));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    A.mu0*=lc; A.mu1*=lc; A.mu2*=lc;
  }
  
  lc = ZZ_p(LeadCoeff(A.nu0));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    A.nu0*=lc; A.nu1*=lc; A.nu2*=lc;
  }
}

// Input: A distinguished infrastructure ideal, A. (rank = 2 only) 
// Output: The ideal A given with the 1-reduced basis.
// See Alg. 4.6 of [LSY] or Alg. 5.1.5 of [L09].

void reduce_basis1(infrastructure_ideal &A){
  int xm, xn, em, en; 
  ZZ_pX temp0, temp1, temp2, k, etamu, etanu;
  ZZ_p lc, lc1;
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);
  ZZ_p _half = to_ZZ_p(half);
  
  //Step 2 in 4.6
  xm = xideg1(A.mu1, A.mu2, A.d);
  xn = xideg1(A.nu1, A.nu2, A.d);

  if(xm < xn){
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0;  A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0 = -temp0; A.nu1 = -temp1; A.nu2 = -temp2;
  }
  if(xm == xn){
    em = etadeg1(A.mu1, A.mu2, A.d); 
    en = etadeg1(A.nu1, A.nu2, A.d);
    if(em < en){
      temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
      A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
      A.nu0 = -temp0; A.nu1 = -temp1; A.nu2 = -temp2;
    }
  }
  
  // Step 3.
  em = etadeg1(A.mu1, A.mu2, A.d); 
  en = etadeg1(A.nu1, A.nu2, A.d);
  if(em >= en){
    // Step 3.1
    while( deg(sqr(A.nu1*_UU*rho) - sqr(A.nu2*_UU2*omega)) - 2*prec  > (double)discdeg(A)/2  ){
      k = divxi1(A);
      temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
      A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
      A.nu0 = k*A.mu0-temp0; A.nu1 = k*A.mu1-temp1; A.nu2 = k*A.mu2-temp2;  
    }
    // Step 3.2
    k = divxi1(A);
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0 = k*A.mu0-temp0; A.nu1 = k*A.mu1-temp1; A.nu2 = k*A.mu2-temp2;

    // Step 3.3
    etamu = A.mu1*_UU*rho - A.mu2*_UU2*omega;
    etanu = A.nu1*_UU*rho - A.nu2*_UU2*omega;
    
    em = deg(etamu); en = deg(etanu);
    if(em==en){
      lc=ZZ_p(LeadCoeff(etamu));
      lc1=ZZ_p(LeadCoeff(etanu));
      lc*=inv(lc1);
      
      A.mu0 -= lc1*A.nu0; 
      A.mu1 -= lc1*A.nu1; 
      A.mu2 -= lc1*A.nu2; 
    }
  }
  //Step 4.
  while(etadeg1(A.mu1, A.mu2, A.d) >= 0){
    k = diveta1(A);
    temp0 = A.nu0; temp1 = A.nu1; temp2 = A.nu2;
    A.nu0 = A.mu0; A.nu1 = A.mu1; A.nu2 = A.mu2;
    A.mu0*=k; A.mu1*=k; A.mu2*=k;
    A.mu2-=temp2; A.mu1-=temp1; A.mu0-=temp0;
  }
  // Step 5.
  if(zetadeg1(A.mu0, A.mu1, A.mu2, A.d) >= 0)
    A.mu0 -= _half*trunczeta1(A.mu0, A.mu1, A.mu2, A.d);
  if(zetadeg1(A.nu0, A.nu1, A.nu2, A.d) >= 0)
    A.nu0 -= _half*trunczeta1(A.nu0, A.nu1, A.nu2, A.d);
  
  // Normalize so mu and nu are monic.
  lc = ZZ_p(LeadCoeff(A.mu0));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    A.mu0*=lc; A.mu1*=lc; A.mu2*=lc;
  }
  
  lc = ZZ_p(LeadCoeff(A.nu0));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    A.nu0*=lc; A.nu1*=lc; A.nu2*=lc;
  }
}

// Input: A distinguished infrastructure ideal, A. (rank = 2 only) 
// Output: The ideal A given with the 2-reduced basis.
// See Alg. 4.6 of [LSY] or Alg. 5.1.5 of [L09].

void reduce_basis2(infrastructure_ideal &A){
  int xm, xn, em, en; 
  ZZ_pX temp0, temp1, temp2, k, etamu, etanu;
  ZZ_p lc, lc1;
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);
  ZZ_p _half = to_ZZ_p(half);

  //Step 2 in 4.6
  xm = xideg2(A.mu1, A.mu2, A.d);
  xn = xideg2(A.nu1, A.nu2, A.d);

  if(xm < xn){
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0;  A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0 = -temp0; A.nu1 = -temp1; A.nu2 = -temp2;
  }
  if(xm == xn){
    em = etadeg2(A.mu1, A.mu2, A.d); 
    en = etadeg2(A.nu1, A.nu2, A.d);
    if(em < en){
      temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
      A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
      A.nu0 = -temp0; A.nu1 = -temp1; A.nu2 = -temp2;
    }
  }

  // Step 3.
  em = etadeg2(A.mu1, A.mu2, A.d); 
  en = etadeg2(A.nu1, A.nu2, A.d);
  if(em >= en){
    // Step 3.1
    while( deg(sqr(A.nu1*_UU2*rho) - sqr(A.nu2*_UU*omega)) - 2*prec  > (double)discdeg(A)/2  ){
      k = divxi2(A);
      temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
      A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
      A.nu0 = k*A.mu0-temp0; A.nu1 = k*A.mu1-temp1; A.nu2 = k*A.mu2-temp2;  
    }
    // Step 3.2
    k = divxi2(A);
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0 = k*A.mu0-temp0; A.nu1 = k*A.mu1-temp1; A.nu2 = k*A.mu2-temp2;

    // Step 3.3
    etamu = A.mu1*_UU2*rho - A.mu2*_UU*omega;
    etanu = A.nu1*_UU2*rho - A.nu2*_UU*omega;
    
    em = deg(etamu); en = deg(etanu);
    if(em==en){
      lc=ZZ_p(LeadCoeff(etamu));
      lc1=ZZ_p(LeadCoeff(etanu));
      lc*=inv(lc1);
      
      A.mu0 -= lc1*A.nu0; A.mu1 -= lc1*A.nu1; A.mu2 -= lc1*A.nu2; 
    }
  }
  //Step 4.
  
  while(etadeg2(A.mu1, A.mu2, A.d) >= 0){
    k = diveta2(A);
    temp0 = A.nu0; temp1 = A.nu1; temp2 = A.nu2;
    A.nu0 = A.mu0; A.nu1 = A.mu1; A.nu2 = A.mu2;
    A.mu0*=k; A.mu1*=k; A.mu2*=k;
    A.mu2-=temp2; A.mu1-=temp1; A.mu0-=temp0;
  }
  // Step 5.
  if(zetadeg2(A.mu0, A.mu1, A.mu2, A.d) >= 0)
    A.mu0 -= _half*trunczeta2(A.mu0, A.mu1, A.mu2, A.d);
  if(zetadeg2(A.nu0, A.nu1, A.nu2, A.d) >= 0)
    A.nu0 -= _half*trunczeta2(A.nu0, A.nu1, A.nu2, A.d);
  
  // Normalize so mu and nu are monic.
  lc = ZZ_p(LeadCoeff(A.mu0));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    A.mu0*=lc; A.mu1*=lc; A.mu2*=lc;
  }
  
  lc = ZZ_p(LeadCoeff(A.nu0));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    A.nu0*=lc; A.nu1*=lc; A.nu2*=lc;
  }
}


// Input: A non-distinguished ideal, A. 
// Output: The ideal A given with the 0-reduced basis.
// See Alg. 4.6 of [LSY] or Alg. 5.1.5 of [L09].

void reduce_basis_nondist(infrastructure_ideal &A){
  int xm, xn, em, en; 
  ZZ_pX temp0, temp1, temp2, k, etamu, etanu;
  ZZ_p lc, lc1;
  cubic_ideal T;
  ZZ_p _half = to_ZZ_p(half);

  xm = xideg(A.mu1, A.mu2, A.d);
  xn = xideg(A.nu1, A.nu2, A.d);
 
  if(xm < xn){
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0;  A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0 = -temp0; A.nu1 = -temp1; A.nu2 = -temp2;
  }
  if(xm == xn){
    em = etadeg(A.mu1, A.mu2, A.d); en = etadeg(A.nu1, A.nu2, A.d);
    if(em < en){
      temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
      A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
      A.nu0 = -temp0; A.nu1 = -temp1; A.nu2 = -temp2;
    }
  }

  // Step 3.
  em = etadeg(A.mu1, A.mu2, A.d); en = etadeg(A.nu1, A.nu2, A.d);
  if(em >= en){
    // Step 3.1
    
    while( deg(sqr(A.nu1*rho) - sqr(A.nu2*omega)) - 2*prec  > (double)discdeg(A)/2  ){
      k = divxi(A);
      temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
      A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
      A.nu0 = k*A.mu0-temp0; A.nu1 = k*A.mu1-temp1; A.nu2 = k*A.mu2-temp2;  
    }
    
    // Step 3.2
    k = divxi(A);
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0 = k*A.mu0-temp0; A.nu1 = k*A.mu1-temp1; A.nu2 = k*A.mu2-temp2;
    
    // Step 3.3
    etamu = A.mu1*rho - A.mu2*omega;
    etanu = A.nu1*rho - A.nu2*omega;
    
    em = deg(etamu); en = deg(etanu);
    if(em==en){
      lc=ZZ_p(LeadCoeff(etamu));
      lc1=ZZ_p(LeadCoeff(etanu));
      lc*=inv(lc1);
      
      A.mu0 -= lc1*A.nu0; A.mu1 -= lc1*A.nu1; A.mu2 -= lc1*A.nu2; 
    }
  } 

  //Step 4.
  while(etadeg(A.nu1, A.nu2, A.d) < 0){
    k = divxi(A);
    
    temp0 = A.mu0; temp1 = A.mu1; temp2 = A.mu2;
    A.mu0 = A.nu0; A.mu1 = A.nu1; A.mu2 = A.nu2;
    A.nu0*=k; A.nu1*=k; A.nu2*=k;
    A.nu2-=temp2; A.nu1-=temp1; A.nu0-=temp0;
  }
  
  while(etadeg(A.mu1, A.mu2, A.d) >= 0){
    k = diveta(A);
    
    temp0 = A.nu0; temp1 = A.nu1; temp2 = A.nu2;
    A.nu0 = A.mu0; A.nu1 = A.mu1; A.nu2 = A.mu2;
    A.mu0*=k; A.mu1*=k; A.mu2*=k;
    A.mu2-=temp2; A.mu1-=temp1; A.mu0-=temp0;
  }

  // Step 5.
  if(zetadeg(A.mu0, A.mu1, A.mu2, A.d) >= 0){
    A.mu0 -= _half*trunczeta(A.mu0, A.mu1, A.mu2, A.d); 
  }
  if(zetadeg(A.nu0, A.nu1, A.nu2, A.d) >= 0){
    A.nu0 -= _half*trunczeta(A.nu0, A.nu1, A.nu2, A.d);
  }

  // Step 6.
  en = etadeg(A.nu1, A.nu2, A.d);
  if((eltdeg(A.nu0, A.nu1, A.nu2, A.d) == en) && (en == 0)){
    A.nu0 -= A.d*LeadCoeff(LeftShift(A.nu0, prec) + A.nu1*rho + A.nu2*omega);
  }

  // Normalize so mu and nu are monic.
  
  lc = ZZ_p(LeadCoeff(A.mu0));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    A.mu0*=lc; A.mu1*=lc; A.mu2*=lc;
  }
  
  lc = ZZ_p(LeadCoeff(A.nu0));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    A.nu0*=lc; A.nu1*=lc; A.nu2*=lc;
  }
}

// Input: A non-distinguished ideal, A. (rank = 1 only) 
// Output: The equivalent distinguished ideal A close to the input.
// Adjusts the distance.
// See Alg. 6.11 of [S01] or Alg. 5.2.4 of [L09].

void reduce_ideal(infrastructure_ideal &A){
  ZZ_pX o0, o1, o2; // temp variables
  ZZ_p i;
  ZZ_pX g;            // gcd of numerator and denominator
  int degtemp = eltdeg(A.mu0, A.mu1, A.mu2, A.d), ed;
  ZZ_pX mu0, mu1, mu2, nu0, nu1, nu2, d;
  ZZ_p _half = to_ZZ_p(half);

  while(degtemp <= 0){
    mu0 = A.mu0, mu1 = A.mu1, mu2 = A.mu2, nu0 = A.nu0, nu1 = A.nu1, nu2 = A.nu2, d = A.d;
    
    A.d0 += degtemp;
  
    // Step 3.4
    // mu = mu^-1 nu = nu*mu^-1
    // Set the denominator.
    A.d = mu1*sqr(mu1)*f;
    A.d += (sqr(mu2*G)*mu2-3*mu0*mu1*mu2*G)*H;
    A.d += mu0*sqr(mu0);
    
    // mu^-1
    o0 = sqr(mu0) - mu1*mu2*G*H;
    o1 = sqr(mu2)*G - mu0*mu1;
    o2 = -mu0*mu2 + sqr(mu1)*H;
    
    A.nu0 = (o1*nu2 + o2*nu1)*G*H + o0*nu0;
    A.nu1 = o0*nu1 + o1*nu0 + o2*nu2*G;
    A.nu2 = o0*nu2 + o1*nu1*H + o2*nu0;
    
    A.mu0 = o0*d; A.mu1 = o1*d; A.mu2 = o2*d;  
    
    // Check to make sure gcd(d, mu, nu) == 1.
    g = GCD(A.d, A.nu2);
    if(g!=1){
      g = GCD(g, A.nu1); g = GCD(g, A.nu0);
      if(g!=1){
	g = GCD(g, A.mu2); g = GCD(g, A.mu1); g = GCD(g, A.mu0);
	if(g!=1) {
	  A.d/=g; A.mu0/=g; A.mu1/=g; A.mu2/=g; A.nu0/=g; A.nu1/=g; A.nu2/=g;
	}
      }
    }

    //Normalize the denominator;
    if(!IsOne(LeadCoeff(A.d))){
      i = inv(LeadCoeff(A.d));
      A.d*=i;
    }

    // Step 3.5 Reduce the ideal {1, mu, nu}.

    reduce_basis_nondist(A);

    degtemp = eltdeg(A.mu0, A.mu1, A.mu2, A.d);
  }

  degtemp = eltdeg(A.nu0, A.nu1, A.nu2, A.d);
  ed = etadeg(A.nu1, A.nu2, A.d);
 
  if( (degtemp < ed) && (ed == 0) ){ 
    mu0 = A.mu0, mu1 = A.mu1, mu2 = A.mu2, nu0 = A.nu0, nu1 = A.nu1, nu2 = A.nu2, d = A.d;
    
    A.d0 += degtemp;
    
    A.d = nu1*sqr(nu1)*f;
    A.d += (sqr(nu2*G)*nu2-3*nu0*nu1*nu2*G)*H;
    A.d += nu0*sqr(nu0);  
   
    //Normalize the denominator;
    if(!IsOne(LeadCoeff(A.d))){
      i = inv(LeadCoeff(A.d));
      A.d*=i;
    } 

    // nu^-1
    o0 = sqr(nu0) - nu1*nu2*G*H;
    o1 = sqr(nu2)*G - nu0*nu1;
    o2 = sqr(nu1)*H - nu0*nu2;
    
    A.mu0 = (o1*mu2 + o2*mu1)*G*H + o0*mu0;
    A.mu1 = o0*mu1 + o1*mu0 + o2*mu2*G;
    A.mu2 = o0*mu2 + o1*mu1*H + o2*mu0;

    // Set (mu, nu) = (mu*nu^{-1}, nu^{-1})
    A.nu0 = o0*d; A.nu1 = o1*d; A.nu2 = o2*d;

    // Check to make sure gcd(d, mu, nu) == 1.
    g = GCD(A.d, A.nu2);
    if(g!=1){
      g = GCD(g, A.nu1); g = GCD(g, A.nu0);
      if(g!=1){
	g = GCD(g, A.mu2); g = GCD(g, A.mu1); g = GCD(g, A.mu0);
	if(g!=1) {
	  A.d/=g; A.mu0/=g; A.mu1/=g; A.mu2/=g; A.nu0/=g; A.nu1/=g; A.nu2/=g;
	}
      }
    }
   
    if(zetadeg(A.nu0, A.nu1, A.nu2, A.d) >= 0){
      A.nu0 -= _half*trunczeta(A.nu0, A.nu1, A.nu2, A.d);      
    }
    
    // Normalize so mu and nu are monic.
    i = ZZ_p(LeadCoeff(A.mu0));
    if(!(IsOne(i)||IsZero(i))){
      i = inv(i);
      A.mu0*=i; A.mu1*=i; A.mu2*=i;
    }
    i = ZZ_p(LeadCoeff(A.nu0));

    if(!(IsOne(i)||IsZero(i))){
      i = inv(i);
      A.nu0*=i; A.nu1*=i; A.nu2*=i;
    }
  }
}

// Input: A non-distinguished ideal, A. (rank = 2 only) 
// Output: The equivalent distinguished ideal A close to the input.
// Adjusts the distance.
// See Alg. 5.2.4 of [L09].

void reduce_ideal_r2(infrastructure_ideal &A){
  ZZ_pX o0, o1, o2; // temp variables
  ZZ_p i;
  ZZ_pX g;            // gcd of numerator and denominator
  int degtemp = eltdeg(A.mu0, A.mu1, A.mu2, A.d), ed;
  ZZ_pX mu0, mu1, mu2, nu0, nu1, nu2, d;
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);
  ZZ_p _half = to_ZZ_p(half);

  while(degtemp <= 0){
    mu0 = A.mu0, mu1 = A.mu1, mu2 = A.mu2, nu0 = A.nu0, nu1 = A.nu1, nu2 = A.nu2, d = A.d;
    
    // Set the distance.
    A.d0 += degtemp;    
    o1 = _UU*A.mu1; o2 = _UU2*A.mu2;
    A.d1 += eltdeg(A.mu0, o1, o2, A.d);
    o1 = _UU2*A.mu1; o2 = _UU*A.mu2;
    A.d2 += eltdeg(A.mu0, o1, o2, A.d);

    // Step 3.4
    // mu = mu^-1 nu = nu*mu^-1
    // Set the denominator.
    A.d = mu1*sqr(mu1)*f;
    A.d += (sqr(mu2*G)*mu2-3*mu0*mu1*mu2*G)*H;
    A.d += mu0*sqr(mu0);
    
    // mu^-1
    o0 = sqr(mu0) - mu1*mu2*G*H;
    o1 = sqr(mu2)*G - mu0*mu1;
    o2 = -mu0*mu2 + sqr(mu1)*H;
    
    A.nu0 = (o1*nu2 + o2*nu1)*G*H + o0*nu0;
    A.nu1 = o0*nu1 + o1*nu0 + o2*nu2*G;
    A.nu2 = o0*nu2 + o1*nu1*H + o2*nu0;
    
    A.mu0 = o0*d; A.mu1 = o1*d; A.mu2 = o2*d;  
    
    // Check to make sure gcd(d, mu, nu) == 1.
    g = GCD(A.d, A.nu2);
    if(g!=1){
      g = GCD(g, A.nu1); g = GCD(g, A.nu0);
      if(g!=1){
	g = GCD(g, A.mu2); g = GCD(g, A.mu1); g = GCD(g, A.mu0);
	if(g!=1) {
	  A.d/=g; A.mu0/=g; A.mu1/=g; A.mu2/=g; A.nu0/=g; A.nu1/=g; A.nu2/=g;
	}
      }
    }

    //Normalize the denominator;
    if(!IsOne(LeadCoeff(A.d))){
      i = inv(LeadCoeff(A.d));
      A.d*=i;
    }

    // Step 3.5 Reduce the ideal {1, mu, nu}.

    reduce_basis_nondist(A);

    degtemp = eltdeg(A.mu0, A.mu1, A.mu2, A.d);
  }

  degtemp = eltdeg(A.nu0, A.nu1, A.nu2, A.d);
  ed = etadeg(A.nu1, A.nu2, A.d);

  if( (degtemp < ed) && (ed == 0) ){ 
    mu0 = A.mu0, mu1 = A.mu1, mu2 = A.mu2, nu0 = A.nu0, nu1 = A.nu1, nu2 = A.nu2, d = A.d;
    
    // Update the distance.
    A.d0 += degtemp;
    o1 = _UU*A.nu1; o2 = _UU2*A.nu2;
    A.d1 += eltdeg(A.nu0, o1, o2, A.d);
    o1 = _UU2*A.nu1; o2 = _UU*A.nu2;
    A.d2 += eltdeg(A.nu0, o1, o2, A.d);
    
    A.d = nu1*sqr(nu1)*f;
    A.d += (sqr(nu2*G)*nu2-3*nu0*nu1*nu2*G)*H;
    A.d += nu0*sqr(nu0);  
    
    //Normalize the denominator;
    if(!IsOne(LeadCoeff(A.d))){
      i = inv(LeadCoeff(A.d));
      A.d*=i;
    } 

    // nu^-1
    o0 = sqr(nu0) - nu1*nu2*G*H;
    o1 = sqr(nu2)*G - nu0*nu1;
    o2 = sqr(nu1)*H - nu0*nu2;
    
    A.mu0 = (o1*mu2 + o2*mu1)*G*H + o0*mu0;
    A.mu1 = o0*mu1 + o1*mu0 + o2*mu2*G;
    A.mu2 = o0*mu2 + o1*mu1*H + o2*mu0;

    // Set (mu, nu) = (mu*nu^{-1}, nu^{-1})
    A.nu0 = o0*d; A.nu1 = o1*d; A.nu2 = o2*d;

    // Check to make sure gcd(d, mu, nu) == 1.
    g = GCD(A.d, A.nu2);
    if(g!=1){
      g = GCD(g, A.nu1); g = GCD(g, A.nu0);
      if(g!=1){
	g = GCD(g, A.mu2); g = GCD(g, A.mu1); g = GCD(g, A.mu0);
	if(g!=1) {
	  A.d/=g; A.mu0/=g; A.mu1/=g; A.mu2/=g; A.nu0/=g; A.nu1/=g; A.nu2/=g;
	}
      }
    }

    if(zetadeg(A.nu0, A.nu1, A.nu2, A.d) >= 0){
      A.nu0 -= _half*trunczeta(A.nu0, A.nu1, A.nu2, A.d);
    }
    
    // Normalize so mu and nu are monic.
    i = ZZ_p(LeadCoeff(A.mu0));
    if(!(IsOne(i)||IsZero(i))){
      i = inv(i);
      A.mu0*=i; A.mu1*=i; A.mu2*=i;
    }
    i = ZZ_p(LeadCoeff(A.nu0));

    if(!(IsOne(i)||IsZero(i))){
      i = inv(i);
      A.nu0*=i; A.nu1*=i; A.nu2*=i;
    }
  }
}

// Input: An infrastructure ideal, A.
// Output: The ideal B, equal to A, given by a canonical basis.
// See Lemma 4.2 of [S01] or Lemma 4.1.2 of [L09].

void inf_to_cubic(infrastructure_ideal &A, cubic_ideal &B){
  ZZ_pX a1, b1, a, b, t, u, w, r1, r2, g;

  B.s = A.d;
  XGCD(B.s2, a1, b1, A.mu2, A.nu2);
  B.s1 = (A.mu1*A.nu2 - A.nu1*A.mu2)/B.s2;

  if(deg(B.s2) == 0){
    t = ZZ_pX();
    a = a1;
    b = b1;
  }
  else{
    XGCD(g, r1, r2, B.s1, B.s2);
    t = r1*(a1*A.mu1 + b1*A.nu1)%B.s2;
    a = a1 - t*A.nu2/B.s2;
    b = b1 + t*A.mu2/B.s2;
  }

  u = B.u = (A.mu0*A.nu2 - A.nu0*A.mu2)/(B.s1*B.s2);
  B.v = (a*A.mu0 + b*A.nu0)/B.s2;
  w = B.w = (a*A.mu1 + b*A.nu1)/B.s2;

  // Normalize s, s', and s".
  B.s*=inv(LeadCoeff(B.s));
  B.s1*=inv(LeadCoeff(B.s1));
  B.s2*=inv(LeadCoeff(B.s2));
  
  B.u %= (B.s/B.s1);
  B.w %= B.s1;
  B.v += u*(B.w-w);
  B.v %= (B.s/B.s2);
}

// Input: An ideal, A, not necessarily primitive
// Output: The primitive ideal B, equivalent to A, given by a canonical basis.
// See Lemma 4.2 of [S01] or Lemma 4.1.2 of [L09].

void inf_to_cubic_gen(infrastructure_ideal &A, cubic_ideal &B){
  ZZ_pX a1, b1, a, b, t, u, w, D;
  ZZ_pX mu0 = A.mu0, mu1 = A.mu1, mu2 = A.mu2, nu0 = A.nu0, nu1 = A.nu1, nu2 = A.nu2, d = A.d;
  
  // Pull out the primitive factor.
  D = GCD(d, nu2);
  if(D!=1){
    D = GCD(D, nu1); D = GCD(D, nu0);
    if(D!=1){
      D = GCD(D, mu2); D = GCD(D, mu1); D = GCD(D, mu0);
      if(D!=1) {
	d/=D; mu0/=D; mu1/=D; mu2/=D; nu0/=D; nu1/=D; nu2/=D;
      }
    }
  }
  
  B.s = d;
  XGCD(B.s2, a1, b1, mu2, nu2);
  B.s1 = (mu1*nu2 - nu1*mu2)/B.s2;

  // Normalize s, s', and s".
  B.s*=inv(LeadCoeff(B.s));
  B.s1*=inv(LeadCoeff(B.s1));
  B.s2*=inv(LeadCoeff(B.s2));

  if( !IsZero(t = (a1*mu1 + b1*nu1)%B.s2) ) {
    if( !IsZero(a = B.s1%B.s2) )
      t /= a;
     else    
       t = ZZ_pX();  // t = 0.
  }
  else
    t = ZZ_pX();  // t = 0.

  a = a1 - t*nu2/B.s2;
  b = b1 + t*mu2/B.s2;
  u = B.u = (mu0*nu2 - nu0*mu2)/(B.s1*B.s2);
  B.v = (a*mu0 + b*nu0)/B.s2;
  w = B.w = (a*mu1 + b*nu1)/B.s2;

  B.u %= (B.s/B.s1);
  B.w %= B.s1;
  B.v += u*(B.w-w);
  B.v %= (B.s/B.s2);

}

// Input: An ideal, A.
// Output: The ideal B, equal to A, given by a (0-)reduced basis.
// See Lemma 4.2 of [S01] or Lemma 4.1.2 of [L09].

void cubic_to_inf(cubic_ideal &A, infrastructure_ideal &B){
  B.d = A.s;
  B.mu0 = A.s1*A.u;
  B.mu1 = A.s1;
  B.mu2 = 0;
  B.nu0 = A.s2*A.v;
  B.nu1 = A.s2*A.w;
  B.nu2 = A.s2;
  reduce_basis_nondist(B);
}

// Input: A distinguished ideal, A.
// Output: The ideal B, equal to A, given by a (0-)reduced basis.
// See Lemma 4.2 of [S01] or Lemma 4.1.2 of [L09].
void cubic_to_inf_dist(cubic_ideal &A, infrastructure_ideal &B){
  B.d = A.s;
  B.mu0 = A.s1*A.u;
  B.mu1 = A.s1;
  B.mu2 = 0;
  B.nu0 = A.s2*A.v;
  B.nu1 = A.s2*A.w;
  B.nu2 = A.s2;
  reduce_basis(B);
}

/* The Baby Step operation in unit rank 1.                          */ 
/* Input: An infrastructure ideal, A, given by a reduced basis.     */
/* Output: The infrastructure ideal, B = bs(A), with reduced basis. */
/* See Algorithm 5.3.14 of [L09].                                   */  

void baby_step_r1(infrastructure_ideal &A, infrastructure_ideal &B){
  ZZ_pX o0, o1, o2; // temp variables
  ZZ_pX g;          // gcd of numerator and denominator

  // Update the degree.
  B.d0 = A.d0 + eltdeg(A.mu0, A.mu1, A.mu2, A.d);

  // Step 3.4
  // mu = mu^-1 nu = nu*mu^-1
  // Set the denominator.
  B.d = power(A.mu0,3)+power(A.mu1,3)*f+power(A.mu2,3)*sqr(G)*H-3*A.mu0*A.mu1*A.mu2*G*H;

  // mu^-1
  o0 = sqr(A.mu0) - A.mu1*A.mu2*G*H;
  o1 = sqr(A.mu2)*G - A.mu0*A.mu1;
  o2 = sqr(A.mu1)*H - A.mu0*A.mu2;

  B.nu0 = (o1*A.nu2 + o2*A.nu1)*G*H + o0*A.nu0;
  B.nu1 = o0*A.nu1 + o1*A.nu0 + o2*A.nu2*G;
  B.nu2 = o0*A.nu2 + o1*A.nu1*H + o2*A.nu0;
    
  B.mu0 = o0*A.d; B.mu1 = o1*A.d; B.mu2 = o2*A.d;  
  
  // Check to make sure gcd(d, mu, nu) == 1.
  g = GCD(B.d, B.nu2);
  if(g!=1){
    g = GCD(g, B.nu1); g = GCD(g, B.nu0);
    if(g!=1){
      g = GCD(g, B.mu2); g = GCD(g, B.mu1); g = GCD(g, B.mu0);
      if(g!=1) {
        B.d/=g; B.mu0/=g; B.mu1/=g; B.mu2/=g; B.nu0/=g; B.nu1/=g; B.nu2/=g;
      }
    }
  }

  //Normalize the basis elements;
  if(!IsOne(LeadCoeff(B.d))){
    B.d *= inv(LeadCoeff(B.d));
  }

  reduce_basis(B);
}


/* The Baby Step operation in the 0-direction in unit rank 2.           */ 
/* Input: An infrastructure ideal, A, given by a 0-reduced basis.       */
/* Output: The infrastructure ideal, B = bs_0(A), with 0-reduced basis. */
/* See Algorithm 5.3.14 of [L09].                                       */  

void baby_step_0_r2(infrastructure_ideal &A, infrastructure_ideal &B){
  ZZ_pX o0, o1, o2; // temp variables
  ZZ_pX a0, a1, a2, t0, t1, t2;
  ZZ_p i;
  ZZ_pX g, autnu;            // gcd of numerator and denominator
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);
  
  // Step 3.4
  // mu = mu^-1 nu = nu*mu^-1
  // Set the denominator.
  if(deg(autnu = LeftShift(A.nu0, prec) + A.nu1*_UU*rho + A.nu2*_UU2*omega) == (prec+deg(A.d))){

    a0 = A.nu0-A.d*ZZ_p(LeadCoeff(autnu)); 
    a1 = A.nu1; a2=A.nu2;

    // Update the distance.
    B.d0 = A.d0 + eltdeg(a0, a1, a2, A.d);
    o1 = _UU*a1; o2 = _UU2*a2;
    B.d1 = A.d1 + eltdeg(a0, o1, o2, A.d);
    o1 = _UU2*a1; o2 = _UU*a2;
    B.d2 = A.d2 + eltdeg(a0, o1, o2, A.d);

    //Step 3.3
    // theta := theta*alpha
    div(o0, (t1*a2 + t2*a1)*G*H + t0*a0, A.d);
    div(o1, t0*a1 + t1*a0 + t2*a2*G, A.d);
    div(o2, t0*a2 + t1*a1*H + t2*a0, A.d);
    t0 = o0; t1 = o1; t2 = o2;
     
    // Step 3.4
    // mu = a^-1 nu = mu*a^-1
    // a^-1
    o0 = sqr(a0) - a1*a2*G*H;
    o1 = sqr(a2)*G - a0*a1;
    o2 = -a0*a2 + sqr(a1)*H;
    
    B.nu0 = (o1*A.mu2 + o2*A.mu1)*G*H + o0*A.mu0;
    B.nu1 = o0*A.mu1 + o1*A.mu0 + o2*A.mu2*G;
    B.nu2 = o0*A.mu2 + o1*A.mu1*H + o2*A.mu0;
    
    B.mu0 = o0*A.d; B.mu1 = o1*A.d; B.mu2 = o2*A.d;
    
    // Set the denominator
    B.d=a1*sqr(a1)*f;
    B.d+=(sqr(a2*G)*a2-3*a0*a1*a2*G)*H;
    B.d+=a0*sqr(a0);
  }
  else{
    // Update the distance.
    B.d0 = A.d0 + eltdeg(A.mu0, A.mu1, A.mu2, A.d);
    o1 = _UU*A.mu1; o2 = _UU2*A.mu2;
    B.d1 = A.d1 + eltdeg(A.mu0, o1, o2, A.d);
    o1 = _UU2*A.mu1; o2 = _UU*A.mu2;
    B.d2 = A.d2 + eltdeg(A.mu0, o1, o2, A.d);

    B.d = power(A.mu0,3)+power(A.mu1,3)*f+power(A.mu2,3)*sqr(G)*H-3*A.mu0*A.mu1*A.mu2*G*H;
    
    // mu^-1
    o0 = sqr(A.mu0) - A.mu1*A.mu2*G*H;
    o1 = sqr(A.mu2)*G - A.mu0*A.mu1;
    o2 = -A.mu0*A.mu2 + sqr(A.mu1)*H;
    
    B.nu0 = (o1*A.nu2 + o2*A.nu1)*G*H + o0*A.nu0;
    B.nu1 = o0*A.nu1 + o1*A.nu0 + o2*A.nu2*G;
    B.nu2 = o0*A.nu2 + o1*A.nu1*H + o2*A.nu0;
    
    B.mu0 = o0*A.d; B.mu1 = o1*A.d; B.mu2 = o2*A.d;  
  }
  
  // Check to make sure gcd(d, mu, nu) == 1.
  g = GCD(B.d, B.nu2);
  if(g!=1){
    g = GCD(g, B.nu1); g = GCD(g, B.nu0);
    if(g!=1){
      g = GCD(g, B.mu2); g = GCD(g, B.mu1); g = GCD(g, B.mu0);
      if(g!=1) {
	B.d/=g; B.mu0/=g; B.mu1/=g; B.mu2/=g; B.nu0/=g; B.nu1/=g; B.nu2/=g;
      }
    }
  }
  
  //Normalize the basis elements;
  if(!IsOne(LeadCoeff(B.d))){
    i = inv(LeadCoeff(B.d));
    B.mu0*=i;  
    B.mu1*=i;
    B.mu2*=i;
    B.nu0*=i;  
    B.nu1*=i;
    B.nu2*=i;
    B.d*=i;
  }
  
  // Step 3.5 Reduce the basis of the ideal {1, mu, nu}.
  reduce_basis(B);
}

/* The Baby Step operation in the 1-direction in unit rank 2.           */ 
/* Input: An infrastructure ideal, A, given by a 1-reduced basis.       */
/* Output: The infrastructure ideal, B = bs_1(A), with 1-reduced basis. */
/* See Algorithm 5.3.14 of [L09].                                       */ 
void baby_step_1_r2(infrastructure_ideal &A, infrastructure_ideal &B){
  ZZ_pX o0, o1, o2; // temp variables
  ZZ_pX a0, a1, a2, t0, t1, t2;
  ZZ_p i;
  ZZ_pX g, autnu; // gcd of numerator and denominator
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);

  // Step 3.4
  // mu = mu^-1 nu = nu*mu^-1
  // Set the denominator.
  if(deg(autnu = LeftShift(A.nu0, prec) + A.nu1*_UU2*rho + A.nu2*_UU*omega) == (prec+deg(A.d))){

    a0 = A.nu0-A.d*ZZ_p(LeadCoeff(autnu));
    a1 = A.nu1; a2=A.nu2;

    // Update the distance.
    B.d0 = A.d0 + eltdeg(a0, a1, a2, A.d);
    o1 = _UU*a1; o2 = _UU2*a2;
    B.d1 = A.d1 + eltdeg(a0, o1, o2, A.d);
    o1 = _UU2*a1; o2 = _UU*a2;
    B.d2 = A.d2 + eltdeg(a0, o1, o2, A.d);

    //Step 3.3
    // theta := theta*alpha
    div(o0, (t1*a2 + t2*a1)*G*H + t0*a0, A.d);
    div(o1, t0*a1 + t1*a0 + t2*a2*G, A.d);
    div(o2, t0*a2 + t1*a1*H + t2*a0, A.d);
    t0 = o0; t1 = o1; t2 = o2;
     
    // Step 3.4
    // mu = a^-1 nu = mu*a^-1
    // a^-1
    o0 = sqr(a0) - a1*a2*G*H;
    o1 = sqr(a2)*G - a0*a1;
    o2 = -a0*a2 + sqr(a1)*H;
    
    B.nu0 = (o1*A.mu2 + o2*A.mu1)*G*H + o0*A.mu0;
    B.nu1 = o0*A.mu1 + o1*A.mu0 + o2*A.mu2*G;
    B.nu2 = o0*A.mu2 + o1*A.mu1*H + o2*A.mu0;
    
    B.mu0 = o0*A.d; B.mu1 = o1*A.d; B.mu2 = o2*A.d;
    
    // Set the denominator
    B.d=a1*sqr(a1)*f;
    B.d+=(sqr(a2*G)*a2-3*a0*a1*a2*G)*H;
    B.d+=a0*sqr(a0);
  }
  else{
    // Update the distance.
    B.d0 = A.d0 + eltdeg(A.mu0, A.mu1, A.mu2, A.d);
    o1 = _UU*A.mu1; o2 = _UU2*A.mu2;
    B.d1 = A.d1 + eltdeg(A.mu0, o1, o2, A.d);
    o1 = _UU2*A.mu1; o2 = _UU*A.mu2;
    B.d2 = A.d2 + eltdeg(A.mu0, o1, o2, A.d);
    
    B.d = power(A.mu0,3)+power(A.mu1,3)*f+power(A.mu2,3)*sqr(G)*H-3*A.mu0*A.mu1*A.mu2*G*H;
    
    // mu^-1
    o0 = sqr(A.mu0) - A.mu1*A.mu2*G*H;
    o1 = sqr(A.mu2)*G - A.mu0*A.mu1;
    o2 = -A.mu0*A.mu2 + sqr(A.mu1)*H;
    
    B.nu0 = (o1*A.nu2 + o2*A.nu1)*G*H + o0*A.nu0;
    B.nu1 = o0*A.nu1 + o1*A.nu0 + o2*A.nu2*G;
    B.nu2 = o0*A.nu2 + o1*A.nu1*H + o2*A.nu0;
    
    B.mu0 = o0*A.d; B.mu1 = o1*A.d; B.mu2 = o2*A.d;  
  }
  
  // Check to make sure gcd(d, mu, nu) == 1.
  g = GCD(B.d, B.nu2);
  if(g!=1){
    g = GCD(g, B.nu1); g = GCD(g, B.nu0);
    if(g!=1){
      g = GCD(g, B.mu2); g = GCD(g, B.mu1); g = GCD(g, B.mu0);
      if(g!=1) {
	B.d/=g; B.mu0/=g; B.mu1/=g; B.mu2/=g; B.nu0/=g; B.nu1/=g; B.nu2/=g;
      }
    }
  }
  
  //Normalize the basis elements;
  if(!IsOne(LeadCoeff(B.d))){
    i = inv(LeadCoeff(B.d));
    B.mu0*=i;  
    B.mu1*=i;
    B.mu2*=i;
    B.nu0*=i;  
    B.nu1*=i;
    B.nu2*=i;
    B.d*=i;
  }
  
  // Step 3.5 1-Reduce the basis of the ideal {1, mu, nu}.
  reduce_basis1(B);
}

/* The Baby Step operation in the 2-direction in unit rank 2.           */ 
/* Input: An infrastructure ideal, A, given by a 2-reduced basis.       */
/* Output: The infrastructure ideal, B = bs_2(A), with 2-reduced basis. */
/* See Algorithm 5.3.14 of [L09].                                       */ 
void baby_step_2_r2(infrastructure_ideal &A, infrastructure_ideal &B){
  ZZ_pX o0, o1, o2; // temp variables
  ZZ_pX a0, a1, a2, t0, t1, t2;
  ZZ_p i;
  ZZ_pX g, autnu;            // gcd of numerator and denominator
  ZZ_p _UU = to_ZZ_p(UU), _UU2 = to_ZZ_p(UU2);

  // Step 3.4
  // mu = mu^-1 nu = nu*mu^-1
  // Set the denominator.
  if(deg(autnu = LeftShift(A.nu0, prec) + A.nu1*rho + A.nu2*omega) == (prec+deg(A.d))){

    a0 = A.nu0-A.d*ZZ_p(LeadCoeff(autnu));
    a1 = A.nu1; a2=A.nu2;

    // Update the distance.
    B.d0 = A.d0 + eltdeg(a0, a1, a2, A.d);
    o1 = _UU*a1; o2 = _UU2*a2;
    B.d1 = A.d1 + eltdeg(a0, o1, o2, A.d);
    o1 = _UU2*a1; o2 = _UU*a2;
    B.d2 = A.d2 + eltdeg(a0, o1, o2, A.d);

    //Step 3.3
    // theta := theta*alpha
    div(o0, (t1*a2 + t2*a1)*G*H + t0*a0, A.d);
    div(o1, t0*a1 + t1*a0 + t2*a2*G, A.d);
    div(o2, t0*a2 + t1*a1*H + t2*a0, A.d);
    t0 = o0; t1 = o1; t2 = o2;
     
    // Step 3.4
    // mu = a^-1 nu = mu*a^-1
    // a^-1
    o0 = sqr(a0) - a1*a2*G*H;
    o1 = sqr(a2)*G - a0*a1;
    o2 = -a0*a2 + sqr(a1)*H;
    
    B.nu0 = (o1*A.mu2 + o2*A.mu1)*G*H + o0*A.mu0;
    B.nu1 = o0*A.mu1 + o1*A.mu0 + o2*A.mu2*G;
    B.nu2 = o0*A.mu2 + o1*A.mu1*H + o2*A.mu0;
    
    B.mu0 = o0*A.d; B.mu1 = o1*A.d; B.mu2 = o2*A.d;
    
    // Set the denominator
    B.d=a1*sqr(a1)*f;
    B.d+=(sqr(a2*G)*a2-3*a0*a1*a2*G)*H;
    B.d+=a0*sqr(a0);
  }
  else{
    B.d = power(A.mu0,3)+power(A.mu1,3)*f+power(A.mu2,3)*sqr(G)*H-3*A.mu0*A.mu1*A.mu2*G*H;
    
    // Update the distance.
    B.d0 = A.d0 + eltdeg(A.mu0, A.mu1, A.mu2, A.d);
    o1 = _UU*A.mu1; o2 = _UU2*A.mu2;
    B.d1 = A.d1 + eltdeg(A.mu0, o1, o2, A.d);
    o1 = _UU2*A.mu1; o2 = _UU*A.mu2;
    B.d2 = A.d2 + eltdeg(A.mu0, o1, o2, A.d);
    
    // mu^-1
    o0 = sqr(A.mu0) - A.mu1*A.mu2*G*H;
    o1 = sqr(A.mu2)*G - A.mu0*A.mu1;
    o2 = -A.mu0*A.mu2 + sqr(A.mu1)*H;
    
    B.nu0 = (o1*A.nu2 + o2*A.nu1)*G*H + o0*A.nu0;
    B.nu1 = o0*A.nu1 + o1*A.nu0 + o2*A.nu2*G;
    B.nu2 = o0*A.nu2 + o1*A.nu1*H + o2*A.nu0;
    
    B.mu0 = o0*A.d; B.mu1 = o1*A.d; B.mu2 = o2*A.d;  
  }
  
  // Check to make sure gcd(d, mu, nu) == 1.
  g = GCD(B.d, B.nu2);
  if(g!=1){
    g = GCD(g, B.nu1); g = GCD(g, B.nu0);
    if(g!=1){
      g = GCD(g, B.mu2); g = GCD(g, B.mu1); g = GCD(g, B.mu0);
      if(g!=1) {
	B.d/=g; B.mu0/=g; B.mu1/=g; B.mu2/=g; B.nu0/=g; B.nu1/=g; B.nu2/=g;
      }
    }
  }
  
  //Normalize the basis elements;
  if(!IsOne(LeadCoeff(B.d))){
    i = inv(LeadCoeff(B.d));
    B.mu0*=i;  
    B.mu1*=i;
    B.mu2*=i;
    B.nu0*=i;  
    B.nu1*=i;
    B.nu2*=i;
    B.d*=i;
  }
  
  // Step 3.5 2-Reduce the basis of the ideal {1, mu, nu}.
  reduce_basis2(B);
}

/* The Giant Step operation in infrastructure.                       */ 
/* Input: Two infrastructure ideals, this, B.                        */
/* Output: The infrastructure ideal, A = this*B, with reduced basis. */
/* See Algorithm 5.3.20 of [L09] or Alg. 7.4 of [S01] for rank = 1.  */ 
infrastructure_ideal infrastructure_ideal::operator* (infrastructure_ideal& B){
  infrastructure_ideal A, result;
  cubic_ideal A2, B2, result2;
  long deg_diff = 0;

  A.mu0 = this->mu0;
  A.mu1 = this->mu1;
  A.mu2 = this->mu2;
  A.nu0 = this->nu0;
  A.nu1 = this->nu1;
  A.nu2 = this->nu2;
  A.d   = this->d;
  A.d0  = this->d0;
  A.d1  = this->d1;
  A.d2  = this->d2;

  inf_to_cubic(A, A2);
  inf_to_cubic(B, B2);

  multiply(result2, A2, B2, deg_diff);
  cubic_to_inf(result2, result);
  
  if(rank == 1){
    result.d0 = A.d0 + B.d0 + deg_diff;
    reduce_basis_nondist(result);
    reduce_ideal(result);
  }
  else{
    result.d0 = A.d0 + B.d0 + deg_diff;
    result.d1 = A.d1 + B.d1 + deg_diff;
    result.d2 = A.d2 + B.d2 + deg_diff;
    reduce_basis_nondist(result);
    reduce_ideal_r2(result);
  }
  return result;
}

/* The Giant Step operation in rank 1 infrastructure.                  */ 
/* Input: Two infrastructure ideals, A and B,                          */
/*          given in reduced and canonical bases, respectively.        */
/*          Bdist - the distance of B as an infrastructure ideal.      */
/* Output: The infrastructure ideal, result = A*B, with reduced basis. */
/* See Algorithm 5.3.20 of [L09] or Alg. 7.4 of [S01].                 */ 
void giant_step_r1(infrastructure_ideal &result, infrastructure_ideal &A, cubic_ideal& B, ZZ &Bdist){
  cubic_ideal A2, result2;
  long deg_diff = 0;

  inf_to_cubic(A, A2);
  
  multiply(result2, A2, B, deg_diff);
  cubic_to_inf(result2, result);
  result.d0 = A.d0 + Bdist + deg_diff;
  reduce_basis_nondist(result);
  reduce_ideal(result);
}

/* The Giant Step operation in rank 2 infrastructure.                     */ 
/* Input: Two infrastructure ideals, A and B,                             */
/*          given in 0-reduced and canonical bases, respectively.         */
/*          Bd0, Bd1, Bd2 - the distance of B as an infrastructure ideal. */
/* Output: The infrastructure ideal, result = A*B, with 0-reduced basis.  */
/* See Algorithm 5.3.20 of [L09].                                         */ 
void giant_step_r2(infrastructure_ideal &result, infrastructure_ideal &A, cubic_ideal& B, ZZ &Bd0, ZZ &Bd1, ZZ &Bd2){
  cubic_ideal A2, result2;
  long deg_diff = 0;

  inf_to_cubic(A, A2);
  
  multiply(result2, A2, B, deg_diff);
  cubic_to_inf(result2, result);
  
  result.d0 = A.d0 + Bd0 + deg_diff;
  result.d1 = A.d1 + Bd1 + deg_diff;
  result.d2 = A.d2 + Bd2 + deg_diff;

  
  reduce_basis_nondist(result);
  reduce_ideal_r2(result);
}

/* The Inverse operation in infrastructure.                              */ 
/* Input: An infrastructure ideals, B, given by a reduced basis.         */
/* Output: The infrastructure ideal, A = Inverse(B), with reduced basis. */
/* See Algorithm 5.3.24 of [L09].                                        */ 
void inverse(infrastructure_ideal &A, infrastructure_ideal &B){
  cubic_ideal A1, B1;
  
  inf_to_cubic(B, B1);
  
  if (IsOne(H)) {  // easy case 
    if (IsOne(B1.s1)) 
      invV1a(A1, B1);
    else 
      invV2a(A1, B1);
  }
  else {
    if (IsOne(B1.s2)) {
      if (IsOne(B1.s1))
	invV1(A1, B1);
      else 
	invV2(A1, B1);
    }
    else 
      genInv(A1, B1);
  }
  
  cubic_to_inf(A1, A);
  if(rank == 1){
    A.d0 = deg(B.d) - B.d0;
    reduce_ideal(A);
  }
  else{
    A.d0 = deg(B.d) - B.d0;
    A.d1 = deg(B.d) - B.d1;
    A.d2 = deg(B.d) - B.d2;
    reduce_ideal_r2(A);
  }
} 

// Regulator extraction in unit rank 1 infrastructure.
// Input: h0 - A multiple of the S-regulator, R^S.
//        bound - a lower bound on R^S.
//        manext - 0 for automatic factoring of h0.
//               - 1 for manual factoring of h0 via an external program.
// Output: The S-regulator, R^S. (returned)
//         factors - the distinct prime factors of h0.
//         exponents - the exponents corresponding to factors.
// See Alg. 6.3.23 of [L09].
ZZ extract(ZZ h0, ZZ bound, vec_ZZ &factors, vec_long &exponents, int manext){
  ZZ h = h0, hx = One;
  ZZ p = One;
  ZZ pe = One;
  ZZ temp;
  int i=0;
  infrastructure_ideal A;

  if(manext){
    cout<<"Input prime factors of h_0, smallest first: [enter '0' when completed]."<<endl;
    
    while(1){
      cin>>p;
      if(IsZero(p))
	break;
      else if(p == h)    
	return hx;
      else{
	factors.SetLength(i+1);
	factors[i] = p;
	i++;
	
      }
    }
  }
  else{
    factor(h, factors, exponents);
  }
  for(i=0; i<factors.length(); i++){
    if(factors[i] < (h/bound)){
      pe = p = factors[i];
      temp = 2*h/pe;
      below(temp, A);
      while(IsOne(A.d) && IsZero(h%pe)){
	pe*=p;
        temp = 2*h/pe;
	below(temp, A);
      }
      hx*=(pe/p);
    }
  }
  return h/hx;
}

// Regulator extraction in unit rank 2 infrastructure.
// Input: h0 - A multiple of the regulator, R.
//        bound - a lower bound on R.
//        manext - 0 for automatic factoring of h0.
//               - 1 for manual factoring of h0 via an external program.
// Output: The regulator, R. (returned)
//         factors - the distinct prime factors of h0.
//         exponents - the exponents corresponding to factors.
// A natural generalization of Alg. 6.3.23 of [L09].
ZZ extract2(ZZ h0, ZZ bound, vec_ZZ &factors, vec_long &exponents, int manext){
  ZZ h = h0, hx = One;
  ZZ p = One;
  ZZ pe = One;
  ZZ temp;
  int i=0;
  infrastructure_ideal A;

  if(manext){
    cout<<"Input prime factors of h_0, smallest first: [enter '0' when completed]."<<endl;
    
    while(1){
      cin>>p;
      if(IsZero(p))
	break;
      else if(p == h)    
	return hx;
      else{
	factors.SetLength(i+1);
	factors[i] = p;
	i++;
	
      }
    }
  }
  else{
    //cout<<"Factoring "<<h<<endl;
    factor(h, factors, exponents);
    //for(i=0; i<factors.length(); i++){
    //  cout<<factors[i]<<" "<<exponents[i]<<"; ";
    //}
    //cout<<endl;
  }
  for(i=0; i<factors.length(); i++){
    if(factors[i] < (h/bound)){
      pe = p = factors[i];
      temp = h/pe;
      //cout<<factors[i]<<" "<<pe<<endl;
      below(temp, Zero, A);
      //cout<<temp<<" "<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
      while(IsOne(A.d) && IsZero(h%(pe*p))){
	pe*=p;
        temp = h/pe;
	//cout<<factors[i]<<" "<<pe<<endl;
	below(temp, Zero, A);
	//cout<<temp<<" "<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
      }
      if(!IsOne(A.d))
	hx*=(pe/p);
      else
	hx*=pe;
    }
  }
  return h/hx;
}

// Integer factorization routine.
// Currently does trial division up to 2111 followed by Pollard's Rho.
// Input: n - an integer
// Output: factors - the distinct prime factors of n
//         exponents - the exponents corresponding to factors.

int factor(ZZ n, vec_ZZ& factors, vec_long& exponents){
  int i;
  int num = 1;
  int exponent = 0;
  ZZ cofactor = n;

  // Extract factors up to 2129 via trial division.
  for(i=0; i<NUMPRIMES; i++){
    while((cofactor%PRIMES[i]) == 0){
      cofactor/=PRIMES[i];
      exponent++;
    }
    if(exponent > 0){
      factors.SetLength(num);
      exponents.SetLength(num);
      factors[num-1] = PRIMES[i];
      exponents[num-1] = exponent;
      num++;
      exponent = 0;
    }
  }
  if(ProbPrime(cofactor)){
    factors.SetLength(num);
    exponents.SetLength(num);
    factors[num-1] = cofactor;
    exponents[num-1] = 1;
    return 1;
  }
  if(IsOne(cofactor))
    return 1;


  // Run Pollard's Rho on the remaining cofactor.
  // Can replace the pollard_rho routine with other factoring algorithms
  // for example: SQUFOF, ECM, etc. if the integers are larger.
  ZZ rhoprime;
  while(!ProbPrime(cofactor)){
    pollard_rho(cofactor, rhoprime);
    cofactor/=rhoprime;
    exponent = 1;
    while((cofactor%rhoprime) == 0){
      cofactor/=rhoprime;
      exponent++;
    }
    factors.SetLength(num);
    exponents.SetLength(num);
    factors[num-1] = rhoprime;
    exponents[num-1] = exponent;
    num++;
    
  }
  factors.SetLength(num);
  exponents.SetLength(num);
  factors[num-1] = cofactor;
  exponents[num-1] = 1;
  return 1;
}

// Pollard's Rho algorithm for factoring.
// A very basic implementation of the algorithm.
// Input: n - an integer
// Output: prime - a prime factor of n.
// See Alg. 8.5.2 of "A Course in Computational Algebraic Number Theory"
//                    by H. Cohen.

void pollard_rho(ZZ n, ZZ& prime){
  ZZ x = to_ZZ(3), y = to_ZZ(3), q = to_ZZ(1), d;
  long i = 1, j = 1; 

  while(1) {
    x  = (x*x - 1)%n;
    y  = (y*y - 1)%n;
    y  = (y*y - 1)%n;
    q *= (x - y); 
    q %= n;
    
    i++;
    if (!j) j=1;
    if ( (i % j) == 0) {
      j++;
      d = GCD(q, n);
      if (!IsOne(d)) {
	prime = d;
	return;
      }
    } // if ( (i % j) == 0)
  } // while n != 1  
}

// Regulator multiple check for unit rank 1 infrastructure.
// Input: h - an integer
// Output: 1 if h is a multiple of R^S.
//         0 if it is not.

int check_inf(ZZ &h){
  infrastructure_ideal D;
  h*=2;
  below(h, D);
  h/=2;
  // This isn't the order. Keep searching.
  if(IsOne(D.d))
    return 1;
  else
    return 0;
}

// Regulator multiple check for unit rank 2 infrastructure.
// Input: h - an integer
// Output: 1 if h is a multiple of R.
//         0 if it is not.

int check_inf2(ZZ &h){
  infrastructure_ideal D;
  //cout<<h<<endl;
  below(h, Zero, D);
  //cout<<D.d0<<" "<<D.d1<<" "<<D.d2<<" "<<D.d<<endl;
  // This isn't the order. Keep searching.
  if(IsOne(D.d))
    return 1;
  else
    return 0;
}
