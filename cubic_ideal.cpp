#ifndef CUBIC_IDEAL_CPP
#define CUBIC_IDEAL_CPP

/* 
   cubic_ideal.cpp

   This file contains functions to perform ideal arithmetic such as 
   multiplication, division, squaring, inversion and reduction. All the 
   special cases have also been implemented with the addition of 
   functions to calculate the minimal element and the canonical basis. 
   For the functions, we will use results and algorithms from:

   [S01] R. Scheidler, "Ideal arithmetic and infrastructure in purely 
                        cubic function fields"
   [B04] M. Bauer, "The Arithmetic of Certain Cubic Function Fields"
   [B05] M. Bauer, Untitled Manuscript
   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)

   There are also miscellaneous functions to calculate the gcd of three 
   elements and bubble sort which were implemented when required by the 
   specific functions. 
   
   How to create the library:
   --------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download cubic_ideal.cpp, cubic_ideal.h and makefile into 
      ../CubicIdeal.
   Compile libcubic via make in ../CubicIdeal.
   
   Make sure that the paths to the NTL libraries are
   specified correctly in the makefile. 
   
   Known Problems: 
   ---------------
   None.
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, ZZ_pEX,  
      GF2X, and GF2EX. 
   2. Include support for characteristic 3 (p=3). See Jonathan Webster.
   3. Include explicit formulas where available: Flon and Oyono, "Fast 
      arithmetic on Jacobians of Picard curves" (deg(G) = 4, H=1)
   4. Include functionality for cubic number fields.
   5. Optimize where possible.

   Authors: 
   Mark Bauer
   Eric Nosal
   Eric Landquist
   Hameeduz Zaman, hameed_uz@hotmail.com, (403) 399-8582
*/

#include "cubic_ideal.h"

NTL_CLIENT



/*******************constructors/destructors**********************/

// A cubic ideal is of the form [s, s'(u + y), s''(v + wp + z)], where
// y = (GH^2)^(1/3) and z = (G^2H)^(1/3).

//default constructor
cubic_ideal::cubic_ideal() {
  s = 1; s1 = 1;  s2 = 1;  u = 0;  v = 0;  w = 0;
}

// constructor
cubic_ideal::cubic_ideal(ZZ_pX& _s, ZZ_pX& _s1, ZZ_pX& _s2, ZZ_pX& _u, ZZ_pX& _v, ZZ_pX& _w) {
  s = _s;  s1 = _s1;  s2 = _s2;  u = _u;  v = _v;  w = _w;
}

// constructor
cubic_ideal::cubic_ideal(const cubic_ideal &A){
  s = A.s;  s1 = A.s1;  s2 = A.s2;  u = A.u;  v = A.v;  w = A.w;
}

//destructor
cubic_ideal::~cubic_ideal() {

}
/*******************constructors/destructors**********************/

/*******************calculator functions************************/

/* 
   This function checks if a given ideal I contains any ramified primes.
   I does not contain any ramified primes if gcd(norm(I), G) == 1.

   -HZ

   I is not ramified iff GCD(s, G*H) = 1.  -- EJL

*/

bool cubic_ideal::is_not_ramified() const {
  if(IsOne(GCD(s,G)) && IsOne(GCD(s, H)))
    return 1; // not ramified
  else  
    return 0; // ramified
}


/* 
   Printing ideals.
   This function prints out the ideal to the screen. Mostly used for debugging
   purposes.

   -HZ
*/

void cubic_ideal::print(){
  cout<<get_s()<<endl;   
  cout<<get_s1()<<endl;  
  cout<<get_s2()<<endl; 
  cout<<get_u()<<endl;   
  cout<<get_v()<<endl; 
  cout<<get_w()<<endl; 
}

/* 
   Ideal validation.
   This function checks if a given ideal is a valid ideal. For more details 
   about the conditions for the validity of an ideal look at [S01, B04, L09].

   -HZ
*/

bool cubic_ideal::is_valid(){
  ZZ_pX temp, r1, r2, modulus;
  ZZ_pX s_G = GCD(get_s(),G);
  ZZ_pX s_H = GCD(get_s(),H);  
  if ((get_s() % (get_s1()*get_s2())) == 0) {    
    if ((H % get_s2()) == 0) {               // do the cheaper statements first
      GCD(temp, get_s()/(s_G * s_H), G*H);  
      if (temp == 1) {
	GCD(temp, get_s1(), H);
	if (temp == 1) {
	  modulus = get_s()/get_s1();
	  rem(r1, H*(get_u()*get_w()-get_v()), modulus); 
	  rem(r2, get_u()*get_u(), modulus);   // reduce function calls here	  
	  if (r1 == r2){
	    modulus = (get_s1()*s_H)/get_s2();
	    ZZ_pX HwP2 = H*get_w()*get_w();
	    rem(r1, HwP2, modulus);
	    rem(r2, get_v(), modulus);
	    if (r1 == r2) {
	      modulus = get_s();
	      rem(r1, H*(G-get_v()*get_w()), modulus);
	      rem(r2, get_u()*(get_v()-HwP2), modulus);
	      if (r1 == r2) {
		return 1;
	      }	  //  if (r1 == r2) {
	    }   //  if (r1 == r2) {
	  }   //  if (r1 == r2){ 
	}   //  if (temp == 1) {
      }   //  if (temp == 1) {
    }   //  if ((get_s2() % get_H()) == 0) {            
  }   //  if ((get_s() % (get_s1()*get_s2())) == 0) {    
  return 0;
}

// For use by MinElt(). Sorts the vector r by weight.

void vecsort(ZZ_pX* r, long* wt){
  long d1, d2, d3; 

  if (IsZero(r[0])) d1=0;
  else d1 = 3*deg(r[0]);
  if (IsZero(r[1])) d2=0;
  else d2 = 3*deg(r[1])+deg(G);
  if (IsZero(r[2])) d3=0;
  else d3 = 3*deg(r[2])+2*deg(G);

  if (d1>d2) {
    if (d1>d3) {
       wt[0]=1;
       wt[1]=d1;
    } else {
        wt[0]=3;
        wt[1]=d3;
    } }
  else {
    if (d2>d3) {
         wt[0]=2;
         wt[1]=d2;
    } else {
     wt[0]=3;
     wt[1]=d3;
     }
  } 
}

// For use by MinElt_general(). Sorts the vector r by weight.

void vecsort_general(ZZ_pX* r, long* wt){
  long d1, d2, d3; 

  if (IsZero(r[0])) d1=0;
  else d1 = 3*deg(r[0]);
  if (IsZero(r[1])) d2=0;
  else d2 = 3*deg(r[1])+deg(f);
  if (IsZero(r[2])) d3=0;
  else d3 = 3*deg(r[2])+2*deg(G) + deg(H);

  if (d1>d2) {
    if (d1>d3) {
       wt[0]=1;
       wt[1]=d1;
    } else {
        wt[0]=3;
        wt[1]=d3;
    } }
  else {
    if (d2>d3) {
         wt[0]=2;
         wt[1]=d2;
    } else {
     wt[0]=3;
     wt[1]=d3;
     }
  } 
}

// Input: An ideal I (H = 1)
// Output: The element, min, of I of minimal norm.
// Algorithm 8.1 of [B04].

void MinElt(ZZ_pX* min, cubic_ideal &I){

  long d1[2],d2[2],d3[2];
  ZZ_pX r, q;
  ZZ_pX r1[3], r2[3], r3[3];

  r1[0] = I.v;
  r1[1] = I.w;
  r1[2] = ZZ_pX(0,1);
  vecsort(r1, d1); 

  r2[0] = I.s1*I.u;
  r2[1] = I.s1;
  r2[2] = ZZ_pX();
  vecsort(r2, d2); 

  r3[0] = I.s;
  r3[1] = ZZ_pX();
  r3[2] = ZZ_pX();
  vecsort(r3, d3); 

  if (d1[1]>d2[1]) {
    swap(r1,r2);
    swap2(d1,d2);
  }
  if (d2[1]>d3[1]) {
    swap(r2,r3);
    swap2(d2,d3);
  }
  if (d1[1]>d2[1]) {
    swap(r1,r2);
    swap2(d1,d2);
  }
   
 while ((d1[0]==d2[0]) || (d1[0]==d3[0]) || (d2[0]==d3[0])) {
   if (d1[0]==d2[0]) {
     DivRem(q,r,r2[d1[0]-1],r1[d1[0]-1]);
     r2[0]=r2[0]-q*r1[0];
     r2[1]=r2[1]-q*r1[1];
     r2[2]=r2[2]-q*r1[2];
     vecsort(r2, d2); 
  }
  else if (d1[0]==d3[0]) {
     DivRem(q,r,r3[d1[0]-1],r1[d1[0]-1]);
     r3[0]=r3[0]-q*r1[0];
     r3[1]=r3[1]-q*r1[1];
     r3[2]=r3[2]-q*r1[2];
     vecsort(r3, d3); 
  }
  else  if (d2[0]==d3[0]) {
     DivRem(q,r,r3[d2[0]-1],r2[d2[0]-1]);
     r3[0]=r3[0]-q*r2[0];
     r3[1]=r3[1]-q*r2[1];
     r3[2]=r3[2]-q*r2[2];
     vecsort(r3, d3); 
  }
  
  if (d1[1]>d2[1]) {
    swap(r1,r2);
    swap2(d1,d2);
  } 

  if (d2[1]>d3[1]) {
    swap(r2,r3);
    swap2(d2,d3);
  }
  
  if (d1[1]>d2[1]) {
    swap(r1,r2);
    swap2(d1,d2);
  }
 }
 min[0] = r1[0];
 min[1] = r1[1];
 min[2] = r1[2];

}

// Input: An ideal I (H != 1)
// Output: The element, min, of I of minimal norm.
// Algorithm 4.5.3 of [L09].

void MinElt_general(ZZ_pX* min, cubic_ideal &I){

  long d1[2],d2[2],d3[2];
  ZZ_pX r, q;
  ZZ_pX r1[3], r2[3], r3[3];

  r1[0] = I.s2*I.v;
  r1[1] = I.s2*I.w;
  r1[2] = I.s2;
  vecsort_general(r1, d1); 

  r2[0] = I.s1*I.u;
  r2[1] = I.s1;
  r2[2] = ZZ_pX();
  vecsort_general(r2, d2); 

  r3[0] = I.s;
  r3[1] = ZZ_pX();
  r3[2] = ZZ_pX();
  vecsort_general(r3, d3); 

  if (d1[1]>d2[1]) {
    swap(r1,r2);
    swap2(d1,d2);
  }
  if (d2[1]>d3[1]) {
    swap(r2,r3);
    swap2(d2,d3);
  }
  if (d1[1]>d2[1]) {
    swap(r1,r2);
    swap2(d1,d2);
  }
  
 while ((d1[0]==d2[0]) || (d1[0]==d3[0]) || (d2[0]==d3[0])) {
   if (d1[0]==d2[0]) {
     DivRem(q,r,r2[d1[0]-1],r1[d1[0]-1]);
     r2[0]=r2[0]-q*r1[0];
     r2[1]=r2[1]-q*r1[1];
     r2[2]=r2[2]-q*r1[2];
     vecsort_general(r2, d2);
  }
  else if (d1[0]==d3[0]) {
     DivRem(q,r,r3[d1[0]-1],r1[d1[0]-1]);
     r3[0]=r3[0]-q*r1[0];
     r3[1]=r3[1]-q*r1[1];
     r3[2]=r3[2]-q*r1[2];
     vecsort_general(r3, d3); 
  }
  else  if (d2[0]==d3[0]) {
     DivRem(q,r,r3[d2[0]-1],r2[d2[0]-1]);
     r3[0]=r3[0]-q*r2[0];
     r3[1]=r3[1]-q*r2[1];
     r3[2]=r3[2]-q*r2[2];
     vecsort_general(r3, d3); 
  }

  if (d1[1]>d2[1]) {
    swap(r1,r2);
    swap2(d1,d2);
  } 

  if (d2[1]>d3[1]) {
    swap(r2,r3);
    swap2(d2,d3);
  }
  
  if (d1[1]>d2[1]) {
    swap(r1,r2);
    swap2(d1,d2);
  }
 }
 min[0] = r1[0];
 min[1] = r1[1];
 min[2] = r1[2];

}


void swap2(long* a, long* b) {
  long temp[2];
  for (int i = 0; i < 2; i++) { temp[i] = a[i]; }
  for (int i = 0; i < 2; i++) { a[i] = b[i]; }
  for (int i = 0; i < 2; i++) { b[i] = temp[i]; }
}

/*
   This function is used by the max_sort function to swap two elements of
   type ZZ_pX. A function to swap long ints already exists in one of the
   NTL libraries (tools.h, I think)

   -HZ
*/

void swap(long* a, long* b) {
  long temp[3];
  for (int i = 0; i < 3; i++) { temp[i] = a[i]; }
  for (int i = 0; i < 3; i++) { a[i] = b[i]; }
  for (int i = 0; i < 3; i++) { b[i] = temp[i]; }
}

/*
   This function is used by the max_sort function to swap two elements of
   type ZZ_pX. A function to swap long ints already exists in one of the
   NTL libraries (tools.h, I think)

   -HZ
*/

void swap(ZZ_pX* a, ZZ_pX* b) {
  ZZ_pX temp[3];
  for (int i = 0; i < 3; i++) { temp[i] = a[i]; }
  for (int i = 0; i < 3; i++) { a[i] = b[i]; }
  for (int i = 0; i < 3; i++) { b[i] = temp[i]; }
}

/* 
   Canonical Basis - Algorithm 0.26 of [B05], Alg. 9.1 of [B04]
   This function takes in the output from MinElt() and 
   finds the canonical basis of the ideal it generates. 
   That is, <min> = I.
   
   -HZ

*/

//  Based on Mark Bauer's code. -- EJL
//  Gives a canonical basis for a principally generated ideal.
//  vector (a,b,c) <-> a+by+cy^2   

void can_basis(cubic_ideal &I, ZZ_pX* min){
  ZZ_pX S, S1, U, V, W, Uprime, Vprime, Wprime;
  ZZ_pX r1, r2, r3, d, q;
  ZZ_pX row1[3], row2[3], row3[3], temp[3];

  row1[0]=min[0];
  row1[1]=min[1];
  row1[2]=min[2];
  
  row2[0]=row1[2]*G; 
  row2[1]=row1[0];
  row2[2]=row1[1];

  row3[0]=row2[2]*G; 
  row3[1]=row2[0];
  row3[2]=row2[1];

  if (row1[2]==0) swap(row1,row2);
  if (row1[2]==0) swap(row1,row3);

  XGCD(d, r1, r2, row1[2], row2[2]);
  temp[2]=d;
  temp[1]=r1*row1[1]+r2*row2[1];
  temp[0]=r1*row1[0]+r2*row2[0];
  r1 = row1[2]/d; r2 = row2[2]/d;
  row2[0]=r1*row2[0]-r2*row1[0];
  row2[1]=r1*row2[1]-r2*row1[1];
  row2[2]=r1*row2[2]-r2*row1[2];
  q=row2[2]/d;
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  row2[0]=row2[0]-q*row1[0];
  row2[1]=row2[1]-q*row1[1];
  row2[2]=row2[2]-q*row1[2];

  XGCD(d, r1, r3, row1[2], row3[2]);
  temp[2]=d;
  temp[1]=r1*row1[1]+r3*row3[1];
  temp[0]=r1*row1[0]+r3*row3[0];
  r1 = row1[2]/d; r3 = row3[2]/d;
  row3[0]=r1*row3[0]-r3*row1[0];
  row3[1]=r1*row3[1]-r3*row1[1];
  row3[2]=r1*row3[2]-r3*row1[2];
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  q=row3[2]/d;
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  row3[0]=row3[0]-q*row1[0];
  row3[1]=row3[1]-q*row1[1];
  row3[2]=row3[2]-q*row1[2];

  if (row2[1]==0) swap(row2,row3);

  XGCD(d, r2, r3, row2[1], row3[1]);
  temp[2]=0;
  temp[1]=d;
  temp[0]=r2*row2[0]+r3*row3[0];
  r2 = row2[1]/d; r3 = row3[1]/d;
  row3[0]=r2*row3[0]-r3*row2[0];
  row3[1]=r2*row3[1]-r3*row2[1];
  //row3[0]=(row2[1]/d)*row3[0]-(row3[1]/d)*row2[0];
  //row3[1]=(row2[1]/d)*row3[1]-(row3[1]/d)*row2[1];
  row3[2]=0;
  row2[0]=temp[0];
  row2[1]=temp[1];
  row2[2]=temp[2];
  q=row3[1]/d;
  row3[0]=row3[0]-q*row2[0];
  row3[1]=0;
  row3[2]=0; 

  // r1 = row1[2]/d; r2 = row2[2]/d
  q=row1[2];
  I.s=row3[0]/q;
  I.s1=row2[1]/q;
  MakeMonic(I.s);
  MakeMonic(I.s1);
  I.u=row2[0]/(I.s1*q);
  I.v=row1[0]/q;
  I.w=row1[1]/q;

  r1=I.s/I.s1;
  rem(I.u, I.u, r1);

  DivRem(q,I.w, I.w, I.s1);  
  I.v=I.v-q*I.u*I.s1;
  rem(I.v, I.v, I.s);

}

/*
   Canonical Basis - Algorithm 0.27 of [B05], Alg. 4.5.4 of [L09]
   H(x) != 1.
   This function takes in the output from MinElt_general() and 
   finds the canonical basis of the ideal it generates. 
   That is, <min> = <D>I.

   -HZ

*/

void can_basis_general(ZZ_pX &D, cubic_ideal &I, ZZ_pX* min){
  ZZ_pX S, S1, U, V, W, Uprime, Vprime, Wprime;
  ZZ_pX r1, r2, r3, d, q;
  ZZ_pX row1[3], row2[3], row3[3], temp[3];

  row1[0]=min[0];
  row1[1]=min[1];
  row1[2]=min[2];
  
  row2[0]=min[2]*G*H; 
  row2[1]=min[0];
  row2[2]=min[1]*H; 

  row3[0]=min[1]*G*H; 
  row3[1]=min[2]*G; 
  row3[2]=min[0];

  if (row1[2]==0) swap(row1,row2);
  if (row1[2]==0) swap(row1,row3);

  XGCD(d, r1, r2, row1[2], row2[2]);
  temp[2]=d;
  temp[1]=r1*row1[1]+r2*row2[1];
  temp[0]=r1*row1[0]+r2*row2[0];
  r1 = row1[2]/d; r2 = row2[2]/d;
  row2[0]=r1*row2[0]-r2*row1[0];
  row2[1]=r1*row2[1]-r2*row1[1];
  row2[2]=r1*row2[2]-r2*row1[2];
  q=row2[2]/d;
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  row2[0]=row2[0]-q*row1[0];
  row2[1]=row2[1]-q*row1[1];
  row2[2]=row2[2]-q*row1[2];

  XGCD(d, r1, r3, row1[2], row3[2]);
  temp[2]=d;
  temp[1]=r1*row1[1]+r3*row3[1];
  temp[0]=r1*row1[0]+r3*row3[0];
  r1 = row1[2]/d; r3 = row3[2]/d;
  row3[0]=r1*row3[0]-r3*row1[0];
  row3[1]=r1*row3[1]-r3*row1[1];
  row3[2]=r1*row3[2]-r3*row1[2];
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  q=row3[2]/d;
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  row3[0]=row3[0]-q*row1[0];
  row3[1]=row3[1]-q*row1[1];
  row3[2]=row3[2]-q*row1[2];

  if (row2[1]==0) swap(row2,row3);

  XGCD(d, r2, r3, row2[1], row3[1]);
  temp[2]=0;
  temp[1]=d;
  temp[0]=r2*row2[0]+r3*row3[0];
  r2 = row2[1]/d; r3 = row3[1]/d;
  row3[0]=r2*row3[0]-r3*row2[0];
  row3[1]=r2*row3[1]-r3*row2[1];
  //row3[0]=(row2[1]/d)*row3[0]-(row3[1]/d)*row2[0];
  //row3[1]=(row2[1]/d)*row3[1]-(row3[1]/d)*row2[1];
  row3[2]=0;
  row2[0]=temp[0];
  row2[1]=temp[1];
  row2[2]=temp[2];
  q=row3[1]/d;
  row3[0]=row3[0]-q*row2[0];
  row3[1]=0;
  row3[2]=0; 

  // r1 = row1[2]/d; r2 = row2[2]/d
  D = GCD(row1[2], row2[1]);

  I.s = row3[0]/D;
  I.s1 = row2[1]/D;
  I.s2 = row1[2]/D;
  MakeMonic(I.s);
  MakeMonic(I.s1);
  MakeMonic(I.s2);
  //MakeMonic(D);
  I.u=row2[0]/row2[1];
  //I.v=row1[0]/row1[2];
  //I.w=row1[1]/row1[2];

  r1=I.s/I.s1;
  rem(I.u, I.u, r1);
  
  XGCD(d, r1, r2, I.s2, I.s/I.s2);
  DivRem(q, I.w, r1*row1[1]/D, I.s1);  
  I.v = r1*row1[0]/D - q*I.u*I.s1;
  rem(I.v, I.v, I.s/I.s2);
}

// Canonical Basis - Algorithm 0.27 of [B05] and Alg. 9.1 of [B04]
// Finds a canonical basis of the ideal, <D>I, generated my the 
// element represented by min.

void can_basis2(ZZ_pX &D, cubic_ideal &I, ZZ_pX* min) {
  ZZ_pX S, S1, U, V, W, Uprime, Vprime, Wprime;
  ZZ_pX r1, r2, r3, d, q;
  ZZ_pX row1[3], row2[3], row3[3], temp[3];

  row1[0]=min[0];
  row1[1]=min[1];
  row1[2]=min[2];
  
  row2[0]=row1[2]*G; 
  row2[1]=row1[0];
  row2[2]=row1[1];

  row3[0]=row2[2]*G; 
  row3[1]=row2[0];
  row3[2]=row2[1];

  row2[0]*=H; 
  row2[2]*=H; 
  row3[2]*=H; 

  if (row1[2]==0) swap(row1,row2);
  if (row1[2]==0) swap(row1,row3);

  XGCD(d, r1, r2, row1[2], row2[2]);
  temp[2]=d;
  temp[1]=r1*row1[1]+r2*row2[1];
  temp[0]=r1*row1[0]+r2*row2[0];
  r1 = row1[2]/d; r2 = row2[2]/d;
  row2[0]=r1*row2[0]-r2*row1[0];
  row2[1]=r1*row2[1]-r2*row1[1];
  row2[2]=r1*row2[2]-r2*row1[2];
  q=row2[2]/d;
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  row2[0]=row2[0]-q*row1[0];
  row2[1]=row2[1]-q*row1[1];
  row2[2]=row2[2]-q*row1[2];

  XGCD(d, r1, r3, row1[2], row3[2]);
  temp[2]=d;
  temp[1]=r1*row1[1]+r3*row3[1];
  temp[0]=r1*row1[0]+r3*row3[0];
  r1 = row1[2]/d; r3 = row3[2]/d;
  row3[0]=r1*row3[0]-r3*row1[0];
  row3[1]=r1*row3[1]-r3*row1[1];
  row3[2]=r1*row3[2]-r3*row1[2];
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  q=row3[2]/d;
  row1[0]=temp[0];
  row1[1]=temp[1];
  row1[2]=temp[2];
  row3[0]=row3[0]-q*row1[0];
  row3[1]=row3[1]-q*row1[1];
  row3[2]=row3[2]-q*row1[2];

  if (row2[1]==0) swap(row2,row3);

  XGCD(d, r2, r3, row2[1], row3[1]);
  temp[2]=0;
  temp[1]=d;
  temp[0]=r2*row2[0]+r3*row3[0];
  r2 = row2[1]/d; r3 = row3[1]/d;
  row3[0]=r2*row3[0]-r3*row2[0];
  row3[1]=r2*row3[1]-r3*row2[1];
  //row3[0]=(row2[1]/d)*row3[0]-(row3[1]/d)*row2[0];
  //row3[1]=(row2[1]/d)*row3[1]-(row3[1]/d)*row2[1];
  row3[2]=0;
  row2[0]=temp[0];
  row2[1]=temp[1];
  row2[2]=temp[2];
  q=row3[1]/d;
  row3[0]=row3[0]-q*row2[0];
  row3[1]=0;
  row3[2]=0; 
  
  // r1 = row1[2]/d; r2 = row2[2]/d
  D=row1[2];
  I.s=row3[0]/D;
  I.s1=row2[1]/D;
  MakeMonic(I.s);
  MakeMonic(I.s1);
  I.u=row2[0]/(I.s1*D);
  I.v=row1[0]/D;
  I.w=row1[1]/D;

  r1=I.s/I.s1;
  rem(I.u, I.u, r1);

  DivRem(q,I.w, I.w, I.s1);  
  I.v=I.v-q*I.u*I.s1;
  rem(I.v, I.v, I.s);
}

/* 
   Ideal reduction Unit Rank 0 - Alg. 10.1 of [B04] (H = 1)
   Algorithm 0.28 of [B05], Alg. 4.5.7 of [L09] (H != 1)
   This function preforms ideal reduction and calls all the necessary 
   functions to do so.

   -HZ
   
   Input:  Nondistinguished ideal I2
   Output: The unique distinguished ideal, I1 equivalent to I2.

   -- EJL
*/

void reduce(cubic_ideal &I1, cubic_ideal &I2){
  ZZ_pX min[3];
  ZZ_pX d;
  cubic_ideal I, J;

  // If this is the identity, then it is obviously reduced.
  if(IsOne(I2.s)){
    I1 = cubic_ideal(I);
    return;
  }

  inverse_nored(I, I2);

  if(IsOne(H)){ 
    MinElt(min, I);
  }
  else{
    MinElt_general(min, I);
  }

  // If the minimum element is 1, then it is already reduced.
  if(IsZero(min[1])&&IsZero(min[2])){
    I1 = cubic_ideal(I2);
    return;
  }

  if(IsOne(H)){ 

    //  if (GCD(min[0], min[1], min[2]) == 1) 
    // then use can basis 1
    if(IsOne(GCD3(min[0], min[1], min[2]))){
      can_basis(J, min);
    }
    
    // else use can basis 2
    else{
      can_basis2(d, J, min);
    }
    
  }
  else{
    can_basis_general(d, J, min);
  }

  if(IsOne(H)){ 

    // I1 = J/I2 -- EJL
    if(J == cubic_ideal()){
      I1 = cubic_ideal(I2);
      return;
    }
    else if(IsOne(d)||IsZero(d)){
      divV3(I1, I, J);
    }
    else{
      divV4(d, I1, I, J);
    }
  }
  else{
    if(J == cubic_ideal()){
      I1 = cubic_ideal(I2);
      return;
    }
    //<D>I = J*I2
    multiply_nored(I1, J, I2);
  }
}

// This function reduces the size of the coefficients of an ideal:
// Old ideal = [s, s'(u  + RHO), s"(v  + w RHO + OMEGA)]
// New ideal = [s, s'(u' + RHO), s"(v' + w'RHO + OMEGA)]
// u' = u (mod s/s'), w' = w (mod s'), and v' = v + u(w' - w) (mod s/s").
//
// Lemma 4.1.6 of [L09], Lemma 4.1.2 of [S01]
// -- EJL

void reduceCoeffs(cubic_ideal &A){
  ZZ_pX u, v, w;

  // Normalize s and s1: Make s and s' monic.
  A.s*=inv(LeadCoeff(A.s));
  A.s1*=inv(LeadCoeff(A.s1));
  A.s2*=inv(LeadCoeff(A.s2));

  u = A.u%(A.s/A.s1);
  w = A.w%A.s1;
  v = (A.v + A.u*(w-A.w))%(A.s/A.s2);

  A.set_u(u);
  A.set_v(v);
  A.set_w(w);

  return;

}

// A kind of constructor.

void ideal(cubic_ideal &I, ZZ_pX S, ZZ_pX S1, ZZ_pX U, ZZ_pX V, ZZ_pX W){
  ZZ_pX vtemp,q, s1u, S2 = ZZ_pX(0,1);
  ZZ_pX redu, redv, redw;

  q = S/S1;
  rem(redu,U,q);

  DivRem(q,redw,W,S1);

  s1u = S1*redu;
  if (IsZero(s1u) || IsZero(S)) vtemp=0;
  else {
    if(deg(S)<=deg(q)) q%=S;
    MulMod(vtemp,s1u,q,S);
  }
  vtemp = V-vtemp;

  rem(redv,vtemp,S); 

  I = cubic_ideal(S, S1, S2, redu, redv, redw);
}


// This generates a random splitting prime cubic_ideal B.
// It assumes that G and H have been initialized. -- EJL
// Algorithm 6.3.18 of [L09], Alg. 0.13 of [B05]

void random(cubic_ideal &A){
  cubic_ideal B, I;
  ZZ_pX P, x, u, d, r1, r2; 
  vec_pair_ZZ_pX_long factors;
  int i, flag=1;

  random(u, 2);
  
  do{
    while(flag){
      //  Generate a random u of degree <2.
      SetCoeff(u, 0, coeff(u,0)+1);
      if(IsZero(coeff(u,0)))
	SetCoeff(u, 1, coeff(u,1)+1);

      // Let s be a factor of u^3-D such that deg(s) > deg(u).
      // Then u^3 = D mod s.
      x = power(u, 3)+f;
      MakeMonic(x);
      berlekamp(factors, x);
      for(i=0; i<factors.length(); i++){
	P = factors[i].a;
	if((deg(P)> deg(u))&&IsOne(GCD(P,H))){
	B.set_s(P);
	  flag=0;
	  break;
	}
      }
    }

    // Reset the flag.
    flag=1;

    B.set_s1(ZZ_pX(0,1));
    B.set_s2(ZZ_pX(0,1));
    B.set_u(u);
    if(deg(H)==0) 
      B.set_v((-sqr(u))%P);
    else{
      XGCD(d, r1, r2, H, P);
      B.set_v((-sqr(u)*r1)%P);
    }
    
    B.set_w(ZZ_pX());

    // If the degree of A <(2/3)g + 1, then it is distinguished/reduced.
    if((deg(f)%3 != 0) && ((double)deg(P) >= ((2.0/3.0)*(double)(deg(G)+deg(H)-1)+1)))
      reduce(A, B);
    else if((deg(f)%3 == 0) && ((double)deg(P) >= ((2.0/3.0)*(double)(deg(G)+deg(H)-2)+1)))
      reduce(A, B);
    else
      A=B;

  } while(A==I);
}


/*
   This function calculates the gcd of three elements, i.e. gcd(a, b, c).
  
   -HZ
*/

ZZ_pX GCD3(const ZZ_pX &A, const ZZ_pX &B, const ZZ_pX &C) {
  ZZ_pX g, g1;
  GCD(g, A, B);
  GCD(g1, g, C);
  return g1;
}


/* 
   Ideal multiplication - Algorithm 0.1 of [B05]
   This function is a specialization of the ideal multiplication function 
   for the case where s1 = 1, s2 = 1.
   
   ATTN: There are three variants (Algs 0.1, 0.2, 0.3) of the algorithm for 
   this case and some testing needs to be done to determine the fastest of 
   the three.
   
   -HZ
*/

void multV1(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX U1, V1, temp, d, r, r1; 

  temp = 0;
  A.set_w(temp);
  
  // half extended Euclidean Algorithm
  // compute r1 such that r1s1 + r2s2 = 1   
  XGCD(d, r1, r, B.s, C.s);
  
  A.set_s(B.s*C.s);  // s = s1*s2

  U1 = (r1*B.s*C.u)+((1-r1*B.s)*B.u);
  rem(temp, U1, A.s);
  A.set_u(temp);
  
  V1 = (r1*B.s*C.v)+((1-r1*B.s)*B.v);
  rem(temp, V1, A.s);
  A.set_v(temp);  
}

/* 
   Ideal multiplication - Algorithm 0.2 of [B05]
   This function is a specialization of the ideal multiplication function 
   for the case where s1 = 1, s2 = 1. 

   ATTN: There are three variants (Algs 0.1, 0.2, 0.3) of the algorithm for 
   this case and some testing needs to be done to determine the fastest of 
   the three.
   
   -HZ
*/

void multV1a(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX U1, V1, temp, d, r, r1; 

  // half extended Euclidean Algorithm
  // compute r1 such that r1s1 + r2s2 = 1   
  XGCD(d, r1, r, B.s, C.s);
  
  A.set_s(B.s*C.s);  // s = s1*s2 

  U1 = (r1*B.s*C.u)+((1-r1*B.s)*B.u);
  rem(temp, U1, A.s);
  A.set_u(temp);
  
  V1 = (B.u*C.u)-(B.u+C.u)*A.u;
  rem(temp, V1, A.s);
  A.set_v(temp);  
}

/*
   Ideal multiplication - Algorithm 0.3 of [B05]
   This function is a specialization of the ideal multiplication function 
   for the case where s1 = 1, s2 = 1. 
   
   ATTN: There are three variants (Algs 0.1, 0.2, 0.3) of the algorithm for 
   this case and some testing needs to be done to determine the fastest of 
   the three.
   
   -HZ
*/

void multV1b(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX U1, V1, temp, d, r, r1; 

  // half extended Euclidean Algorithm
  // compute r1 such that r1s1 + r2s2 = 1   
  XGCD(d, r1, r, B.s, C.s);
  
  A.set_s(B.s*C.s);  // s = s1*s2

  U1 = (r1*B.s*C.u)+((1-r1*B.s)*B.u);
  rem(temp, U1, A.s);
  A.set_u(temp);
  
  V1 = A.u*A.u; 
  rem(temp, V1, A.s);
  A.set_v(temp); 
}


/*
   Ideal multiplication - Algorithm 0.4 of [B05]
   This function is a specialization of the ideal multiplication function 
   for the case where s1 != 1, s2 = 1. 
   
   -HZ
*/

void multV2(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX U1, W1, V1, temp, d, r, r1, r2, q;

  // half extended Euclidean Algorithm
  // compute r1 such that r1s1 + r2s2 = 1   
  XGCD(d, r1, r, B.s, C.s);
  
  rem(r2, r1*B.s1, C.s/C.s1);
  
  // set S
  A.set_s(B.s*C.s);
  A.set_s1(B.s1*C.s1);
  
  U1 = (r2*(B.s/B.s1)*C.u)+(1-r2*(B.s/B.s1))*B.u;
  rem(temp, U1, A.s/A.s1);
  A.set_u(temp);
  
  W1 = (r1*B.s*C.w)+(1-r1*B.s)*B.w;
  V1 = (r1*B.s*C.v)+(1-r1*B.s)*B.v;
  
  div(q, W1, A.s1);
  rem(temp, W1, A.s1);
  A.set_w(temp);
  
  rem(temp, V1-q*A.s1*A.u, A.s);
  A.set_v(temp);  
}

/*
   Ideal multiplication - Algorithm 0.5 of [B05], Alg. 4.4.7 of [L09]
   This function performs ideal multiplication for the general case. 

   -HZ
*/

void genMult(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX U1, W1, V1, temp, d, r, r1, r2, q; 

  // half extended Euclidean Algorithm
  // compute r1 such that r1s1 + r2s2 = 1   
  XGCD(d, r1, r, B.s, C.s);
  
  rem(r2, r1*B.s1, C.s/C.s1);
  
  // set S
  A.set_s(B.s*C.s);
  A.set_s1(B.s1*C.s1);
  A.set_s2(B.s2*C.s2);
  
  U1 = (r2*(B.s/B.s1)*C.u)+(1-r2*(B.s/B.s1))*B.u;
  rem(temp, U1, A.s/A.s1);
  A.set_u(temp);
  
  rem(r2, r1*B.s2, C.s/C.s2);

  V1 = (r2*(B.s/B.s2)*C.v)+(1-r2*(B.s/B.s2))*B.v;
  W1 = (r2*(B.s/B.s2)*C.w)+(1-r2*(B.s/B.s2))*B.w;
  
  div(q, W1, A.s1);
  rem(temp, W1, A.s1);
  A.set_w(temp);
  
  rem(temp, V1-q*A.s1*A.u, A.s/A.s2);
  A.set_v(temp);  
}

/* 
   Ideal multiplication - Algorithm 4.4.9 of [L09]
   This function performs ideal multiplication for the case where the
   product is primitive.
   A = B*C, H(x) != 1.

   -HZ
   Modified by EJL
*/

void mult29(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX U1, U2, W1, V1, temp, r, r1, r2, r3, r4, q; 
  ZZ_pX D, d1, dH, d2;
  ZZ_pX g, g1;
  
  XGCD(D, r1, r2, B.s/B.s1, C.s/C.s1);
  d1 = GCD(D, B.u-C.u)/GCD(D, G*H); 
  GCD(dH, D, H);

  // set S
  A.set_s(B.s*C.s*d1/D);
  A.set_s1(B.s1*C.s1*D/(d1*dH));
  A.set_s2(B.s2*C.s2*dH);
  
  U1 = B.u+((C.u-B.u)*(r1*B.s/(B.s1*D)));
  XGCD(d2, r3, r4, 3*sqr(U1), d1);
  U2 = U1 - r3*((power(U1,3) + f)/d2);
  rem(temp, U2, A.s/A.s1); 
  A.set_u(temp);
  
  // six element GCD  
  XGCD(g, r1, r2, B.s*C.s2,         B.s2*C.s);
  V1 =         r1*B.s*C.s2*C.v + r2*B.s2*C.s*B.v;
  W1 =         r1*B.s*C.s2*C.w + r2*B.s2*C.s*B.w;
  
  if(g != A.s2){
    XGCD(g1, r1, r2, g, B.s1*C.s1*H);
    V1 = r1*V1 +     r2*B.s1*C.s1*B.u*C.u;
    W1 = r1*W1 +     r2*B.s1*C.s1*(B.u+C.u);
    g = g1;
  
    if(g != A.s2){
      XGCD(g1, r1, r2, g, B.s1*C.s2*(B.u+C.w*H)); 
      V1 = r1*V1 +     r2*B.s1*C.s2*(B.u*C.v+G*H); 
      W1 = r1*W1 +     r2*B.s1*C.s2*(B.u*C.w+C.v);
      g=g1;
      
      if(g != A.s2){
	XGCD(g1, r1, r2, g, C.s1*B.s2*(C.u+B.w*H));  
	V1 = r1*V1 +     r2*C.s1*B.s2*(C.u*B.v+G*H); 
	W1 = r1*W1 +     r2*C.s1*B.s2*(C.u*B.w+B.v);
	g=g1;
      
	if(g != A.s2){
	  XGCD(g1, r1, r2, g, B.s2*C.s2*(B.v+C.v+B.w*C.w*H)); 
	  V1 = r1*V1 +     r2*B.s2*C.s2*(B.v*C.v+(B.w+C.w)*G*H); 
	  W1 = r1*W1 +     r2*B.s2*C.s2*(B.w*C.v+B.v*C.w+G); 
	}
      }
    }
  }
  
  if(IsOne(A.s2)){
    rem(A.w, W1, A.s1);
    div(q, A.w-W1, A.s1);
    A.set_v((V1 + q*A.s1*A.u)%A.s);
  }
  else{
    V1/=A.s2;
    if(IsZero(W1%A.s2)){    
      W1/=A.s2;
      
      A.w = W1%A.s1;
      A.v = (V1+A.u*(A.w - W1))%(A.s/A.s2);
    }
    else{
      if(IsOne(A.s1)){
	A.v = (V1 - (A.u/A.s2)*W1)%(A.s/A.s2);
	A.w = ZZ_pX();
      }
      else{
	XGCD(g1, r1, r2, A.s2, A.s1);
	A.w = W1*r1%A.s1;
	A.v = (V1 + (A.u/A.s2)*(A.s2*A.w - W1))%(A.s/A.s2);
      }
    }
  }

}

/* 
   Ideal multiplication - Algorithm 0.30 of [B05], Alg. 7.2 of [B04]
   This function performs ideal multiplication for the case where the
   product is primitive, [Ba] case. (H = 1)

   -HZ

   All fixed.
   A = B*C, H(x) = 1 -- EJL

*/

void mult30(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX u, W1, V1, temp, r, r1, r2, q; 
  ZZ_pX d, d1, d2;
  ZZ_pX g, g1;
  
  // Step 1
  XGCD(d, r1, r2, B.s/B.s1, C.s/C.s1);
  
  // Step 2
  d1 = GCD(d, B.u-C.u)/GCD(d, G); 
  
  // Step 3
  A.set_s(B.s*C.s*d1/d);
  A.set_s1(B.s1*C.s1*d/d1);
  u = B.u+((C.u-B.u)*(r1*B.s/(B.s1*d)));

  // Step 4
  if(!IsOne(d1)){
    XGCD(d2, r1, r2, 3*sqr(u), d1);
    
    // Step 5 and 6
    u = u - r1*((u*sqr(u)+G)/d2); 
    rem(A.u, u, A.s/A.s1); 
  }
  else
    rem(A.u, u, A.s/A.s1);

  // Step 7
  // six element GCD  
  XGCD(g, r1, r2, B.s, B.v+C.v+B.w*C.w);

  V1 = r1*B.s*C.v + r2*(B.v*C.v+(B.w+C.w)*G); 
  W1 = r1*B.s*C.w + r2*(B.w*C.v+B.v*C.w+G); 

  if(!IsOne(g)){
    XGCD(g1, r1, r2, g ,B.s1*C.s1); 
    V1 = r1*V1 + r2*B.s1*C.s1*B.u*C.u;
    W1 = r1*W1 + r2*B.s1*C.s1*(B.u+C.u);
    g = g1;

    if(!IsOne(g)){
      XGCD(g1, r1, r2, g , B.s1*(B.u+C.w));
      V1 = r1*V1 + r2*B.s1*(B.u*C.v+G); 
      W1 = r1*W1 + r2*B.s1*(B.u*C.w+C.v);
      g=g1;
    
      if(!IsOne(g)){
	XGCD(g1, r1, r2, g , C.s1*(C.u+B.w));
	V1 = r1*V1 + r2*C.s1*(C.u*B.v+G); 
	W1 = r1*W1 + r2*C.s1*(C.u*B.w+B.v);
	g=g1;
      
	if(!IsOne(g)){
	  XGCD(g1, r1, r2, g , C.s);
	  V1 = r1*V1 + r2*C.s*B.v;
	  W1 = r1*W1 + r2*C.s*B.w;
	  g=g1;
	}
      }
    }
  }

  // Step 10.
  rem(A.w, W1, A.s1);
  div(q, A.w-W1, A.s1);
  
  // Step 11
  rem(A.v, V1 + q*A.s1*A.u, A.s);

}

// Algorithm 0.30 of [B05]
// Simplified for the case that B.s1 = C.s1 = 1.
// A = B*C

void mult30s(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C, ZZ_pX &DD, ZZ_pX &r){

  ZZ_pX D1, DG, r1, r2, q;
  ZZ_pX temp;
  ZZ_pX Uprime, Vprime, Wprime;
  ZZ_pX one = ZZ_pX(0,1);
  cubic_ideal X;

  GCD(D1, DD, B.u-C.u);

  if (!IsOne(D1)){
    GCD(DG, DD, G); 
    D1=D1/DG;
  }

  X.s  = B.s*(C.s/DD)*D1;
  X.s1 = DD/D1;

  Uprime = B.u-(B.u-C.u)*r*(B.s/DD);

  if (!IsOne(D1)) {  
    XGCD(X.u, r1, r2, D1, 3*sqr(Uprime)); 
    Uprime = Uprime - r2*(power(Uprime,3) + G)/X.u; 
  }

  rem(X.u, Uprime, X.s/X.s1);

  Vprime = B.u*C.u;
  Wprime = B.u+C.u;

  div(q, Wprime, X.s1);
  rem(X.w, Wprime, X.s1);
  
  rem(X.v, Vprime - q*X.s1*X.u, X.s);
}
/* 
   Ideal multiplication - Algorithm 4.4.11  of [L09] - Modified via EJL
   This function performs ideal multiplication for the case where the
   product is non-primitive.

   -HZ
   
   <D>A = B*C,  H(x) != 1.
   -- EJL
   

*/

void mult31(ZZ_pX &D,ZZ_pX &DH,ZZ_pX &D1,ZZ_pX &D2,ZZ_pX &D3,ZZ_pX &D4,ZZ_pX &D5, cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  cubic_ideal I1, I2, I3, I4;
  
  I1.set_s(B.s/D);
  I1.set_s1(B.s1/(D2*D3));
  I1.set_s2(B.s2/(D5*DH));
  I1.set_u(B.u);
  I1.set_v(B.v);
  I1.set_w(B.w);
  reduceCoeffs(I1);

  I2.set_s(C.s/D);
  I2.set_s1(C.s1/(D1*D3));
  I2.set_s2(C.s2/(D4*DH));
  I2.set_u(C.u);
  I2.set_v(C.v);
  I2.set_w(C.w);
  reduceCoeffs(I2);

  I3.set_s(D3*DH);
  I3.set_u((B.w+C.w)*H); 
  I3.set_v(-B.w*C.w*H); 

  if(!IsOne(I3.s)){
    reduceCoeffs(I3);
    mult29(I4, I1, I3);
    mult29(A, I4, I2); 
  }
  else
    mult29(A, I1, I2); 
}

/* 
   Ideal multiplication - Algorithm 0.32 of [B05] (Lemma 7.3 of [B04])
   This function performs ideal multiplication for the case where the
   product is non-primitive, [Ba] case. H = 1.

  -HZ

  <D>A = B*C -- EJL

*/


void mult32(ZZ_pX &D, ZZ_pX &D1, ZZ_pX &D2, ZZ_pX &D3, cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  cubic_ideal I1, I2, I3, I4;

  I1.set_s(B.s/D);
  I1.set_s1(B.s1/(D2*D3));
  I1.set_u(B.u);
  I1.set_v(B.v);
  I1.set_w(B.w);

  reduceCoeffs(I1);

  I2.set_s(C.s/D);
  I2.set_s1(C.s1/(D1*D3));
  I2.set_u(C.u);
  I2.set_v(C.v);
  I2.set_w(C.w);

  reduceCoeffs(I2);
 
  if(IsOne(D3))
    mult30(A, I1, I2);
  else{
    I3.set_s(D3);
    I3.set_u(B.w+C.w);
    I3.set_v(-B.w*C.w);
    reduceCoeffs(I3);
    
    mult30(I4, I1, I2); 
    mult30(A, I3, I4); 
  }
}

/* 
   Ideal squaring 
   This function selects the appropriate specialization of the square 
   function to call. (adapted from existing code from Dr. Bauer)  
   
   -HZ

   Finds A ~ B^2 -- EJL

*/

void square(cubic_ideal &A, const cubic_ideal &B) {
  cubic_ideal X;

  if(B.is_not_ramified()) { // No ramified primes
    if(IsOne(B.s1)){
	sqrV1(X, B);
    }
    else {
      sqrV2(X, B);
    }
  }
  else {
    genSqr(X, B);
  }
  reduce(A, X);
  return;
}

// Squaring for infrastructure. <D>A = B^2

void square(cubic_ideal &A, const cubic_ideal &B, long &ddiff){ 
  if(B.is_not_ramified()) { // No ramified primes

    ddiff = 0;
    if(IsOne(B.s1)){
	sqrV1(A, B);
    }
    else {
      sqrV2(A, B);
    }
  }
  else {
    ddiff = deg(GCD(B.s1, G)) + deg(B.s2) - deg(GCD(B.s, G*H));
    genSqr(A, B);
  }
}

/*
   Ideal squaring - Algorithm 0.7 of [B05]
   This function is a specialization of the ideal squaring function 
   for the case where s1 = s2 = 1, no ramified primes. 
   
   -HZ
*/

void sqrV1(cubic_ideal &A, const cubic_ideal &B) {
  ZZ_pX temp, d, r1, r2;

  A.set_s(sqr(B.s));
 
  if(IsOne(H)){
    XGCD(d, r1, r2, 3*sqr(B.u), B.s);
    A.set_u((B.u - r1*(power(B.u,3)+G))%A.s); 
    A.set_v((sqr(B.u)-2*B.u*A.u)%A.s);
  }
  else{
    XGCD(d, r1, r2, 3*sqr(B.u), B.s);
    A.set_u((B.u - r1*(power(B.u,3) + f))%A.s);
    
    XGCD(d, r1, r2, H, A.s); 
    A.set_v(((sqr(B.u)-2*B.u*A.u)*r1)%A.s);
    
  }
}  


/*
   Ideal squaring - Algorithm 0.8 of [B05]
   This function is a specialization of the ideal squaring function 
   for the case where s1 = 1, s2 != 1, no ramified primes. 
   
   ATTN: Need to test whether Algorithm 0.8 or 0.9 is faster.
   
   -HZ
*/

void sqrV1a(cubic_ideal &A, const cubic_ideal &B) {
  ZZ_pX temp, d, r1, r2;

  A.set_s(sqr(B.s));
 
  XGCD(d, r1, r2, 3*sqr(B.u), B.s);
  A.set_u((B.u - r1*(power(B.u,3)+G))%A.s); 
  
  A.set_v((sqr(B.u)-2*B.u*A.u)%A.s);
}

/*
   Ideal squaring - Algorithm 0.9 of [B05]
   This function is a specialization of the ideal squaring function 
   for the case where s1 = 1, s2 = 1, no ramified primes. 
   
   ATTN: Need to test whether Algorithm 0.8 or 0.9 is faster.
   
   -HZ
*/

void sqrV1b(cubic_ideal &A, const cubic_ideal &B) {
  ZZ_pX temp, d, r1, r2;

  A.set_s(sqr(B.s));
 
  XGCD(d, r1, r2, 3*sqr(B.u), B.s);
  A.set_u((B.u - r1*(power(B.u,3)+G))%A.s); 
  
  XGCD(d, r1, r2, H, A.s); 
  A.set_v((-1*A.u*A.u*r1)%A.s);

}


/*
   Ideal squaring - Algorithm 0.10 of [B05]
   This function is a specialization of the ideal squaring function 
   for the case where s1 != 1, s2 = 1, no ramified primes.
  
   -HZ
*/

void sqrV2(cubic_ideal &A, const cubic_ideal &B) {
  ZZ_pX V1, temp, d, z, r1, r2, wcubed, tv = B.v, tw = B.w;

  A.set_s(sqr(B.s));
  A.set_s1(sqr(B.s1));

  if(IsOne(H)){
    temp = 2*tv+sqr(tw);
    
    XGCD(d, r1, z, B.s, temp);
    while(!IsOne(d)){
      tw += B.s1;
      tv += B.u*B.s1;
      temp = 2*tv+sqr(tw);
      XGCD(d, r1, z, B.s, temp);
    }

    wcubed= power(tw,3);

    XGCD(d, r1, r2, 3*sqr(B.u), B.s/B.s1);
    A.set_u((B.u - r1*(power(B.u,3)+G))% (A.s/A.s1) );
    A.set_w((tw - z*(wcubed - G))%A.s1 );
    V1 = tv + A.u*(A.w-tw) + z*(A.u*(wcubed - G) + 2*G*tw-tv*(tv+sqr(tw)));
    A.set_v(V1 % A.s);
  }
  else {
    temp = 2*tv+H*sqr(tw);     
    XGCD(d, r1, z, B.s, temp);
    while(!IsOne(d)){
      tw += B.s1;
      tv += B.u*B.s1;
      temp = 2*tv+H*sqr(tw); 
      XGCD(d, r1, z, B.s, temp);
    }
    wcubed= power(tw,3);

    XGCD(d, r1, r2, 3*sqr(B.u), B.s/B.s1);
    A.set_u((B.u - r1*(power(B.u,3) + f))% (A.s/A.s1) );
    
    A.set_w((tw - z*(H*wcubed - G))%A.s1 ); 
    V1 = tv + A.u*(A.w-tw) + z*(A.u*(H*wcubed - G) + 2*G*H*tw-tv*(tv+H*sqr(tw)));
    A.set_v(V1 % A.s);
  }
}


/*
   Ideal squaring - Algorithm 0.11 of [B05]
   Algorithm 4.4.8 of [L09]
   This function performs ideal squaring for the general case. 
   
   -HZ
*/

void genSqr(cubic_ideal &A, const cubic_ideal &B) {
  ZZ_pX U1, W1, V1, temp, temp2, d, y, z, r1, r2, sa, sb, wcubed, tv = B.v, tw = B.w;
  ZZ_pX s_G = GCD(B.s,G); 
  ZZ_pX s1_G = GCD(B.s1,G); 
  ZZ_pX s_H = GCD(B.s,H);  
  ZZ_pX sonsgh = B.s/(s_G*s_H);

  if(IsOne(H)){
    temp2 = 2*tv+sqr(tw);
    
    XGCD(d, r1, z, sonsgh, temp2);
    while(!IsOne(d)){
      tw += B.s1;
      tv += B.u*B.s1;
      temp2 = 2*tv+sqr(tw);
      XGCD(d, r1, z, sonsgh, temp2);
    }
    wcubed = power(tw,3);

    A.set_s(B.s*sonsgh);
    A.set_s1((sqr(B.s1)*s_G)/power(s1_G,3));
    A.set_s2(s_H/B.s2);
    
    XGCD(d, y, r2, 3*sqr(B.u), B.s*s1_G/(s_G*s_H*B.s1));
    XGCD(d, r1, r2, sa = s_H*s1_G, sb = sqr(B.s*s1_G/(s_G*s_H*B.s1)));
    U1 = (B.u - y*(power(B.u, 3)+G))*r1*sa % (sa*sb); 
    rem(temp, U1, A.s/A.s1);
    A.set_u(U1);
    
    XGCD(d, r1, r2, sa = s_G/s1_G, sb = sqr(B.s1/s1_G));
    W1 = (tw - z*(wcubed - G))*r1*sa % (sa*sb); 
    A.set_w(W1%A.s1);
    
    V1 = tv + A.u*(A.w-tw) + z*(A.u*(wcubed - G) + 2*G*tw-tv*(tv+sqr(tw)));
    XGCD(d, r1, r2, sa = s_G*B.s2, sb = sqr(sonsgh));
    rem(temp, V1*r1*s_G*B.s2 % (sa*sb), A.s);
    A.set_v(temp);
  }
  else{
    temp2 = 2*tv+H*sqr(tw);
    XGCD(d, r1, z, sonsgh, temp2);
    while(!IsOne(d)){
      tw += B.s1;
      tv += B.u*B.s1;
      temp2 = 2*tv+H*sqr(tw); 
      XGCD(d, r1, z, sonsgh, temp2);
    }
    wcubed = power(tw,3);
    
    A.set_s(B.s*sonsgh);
    A.set_s1((sqr(B.s1)*s_G)/power(s1_G,3));
    A.set_s2(s_H/B.s2);
    
    XGCD(d, y, r2, 3*sqr(B.u), B.s*s1_G/(s_G*s_H*B.s1));
    XGCD(d, r1, r2, sa = s_H*s1_G, sb = sqr(B.s*s1_G/(s_G*s_H*B.s1)));
    U1 = (B.u - y*(power(B.u, 3) + f))*r1*sa % (sa*sb);
    rem(temp, U1, A.s/A.s1);
    A.set_u(U1);
    
    XGCD(d, r1, r2, sa = s_G/s1_G, sb = sqr(B.s1/s1_G));
    W1 = (tw - z*(H*wcubed - G))*r1*sa % (sa*sb);
    A.set_w(W1%A.s1);
    
    V1 = tv + A.u*(A.w-tw) + z*(A.u*(H*wcubed - G) + 2*G*H*tw-tv*(tv+H*sqr(tw)));
    XGCD(d, r1, r2, sa = s_G*B.s2, sb = sqr(sonsgh));
    rem(temp, V1*r1*s_G*B.s2 % (sa*sb), A.s);
    A.set_v(temp);
  }
}

/*
   Ideal inversion 
   This function selects the appropriate specialization of the inverse 
   function to call. (adapted from existing code from Dr. Bauer)  
   
   -HZ

   A = <d>B^{-1}, for some d, where A is the primitive distinguished
   ideal in the class of B^{-1}. -- EJL

*/

void inverse(cubic_ideal &A, cubic_ideal &B) {
  cubic_ideal Inv, temp;

  if (IsOne(H)) {  // easy case if (IsOne(B.H)) {
    if (IsOne(B.s1))
      invV1a(temp, B);
    else
      invV2a(temp, B);
  }
  else {
    if (IsOne(B.s2)) {
      if (IsOne(B.s1))
	invV1(temp, B);
      else 
	invV2(temp, B);
    }
    else
      genInv(temp, B);
  }  
  reduce(A, temp);
}

/*
   Ideal inversion - without reduction: for the reduction function.
   This function selects the appropriate specialization of the inverse 
   function to call. (adapted from existing code from Dr. Bauer)  
   
   -HZ

   A = <d>B^{-1}, for some d, where A is the primitive distinguished
   ideal in the class of B^{-1}. -- EJL

*/

void inverse_nored(cubic_ideal &A, cubic_ideal &B) {
  if (IsOne(H)) {  // easy case 
    if (IsOne(B.s1)) 
      return invV1a(A, B);
    else 
      return invV2a(A, B);
  }
  else {
    if (IsOne(B.s2)) {
      if (IsOne(B.s1)) 
	return invV1(A, B);
      else 
	return invV2(A, B);
    }
    else 
      return genInv(A, B);
  }
}
/*
   Ideal inversion - Algorithm 0.16 of [B05]
   This function is a specialization of the ideal inversion function 
   for the case where s1 = 1, s2 = 1. H != 1.
   A = <s>B^{-1}

   -HZ
*/

void invV1(cubic_ideal &A, cubic_ideal &B) {
  ZZ_pX d, r1, r2;
  ZZ_pX s_H = GCD(B.s, H); 

  A.set_s(B.s);
  A.set_s1(B.s/s_H);  
  A.set_s2(s_H);
 
  XGCD(d, r1, r2, H, A.s1); 
  
  A.set_w((-B.u*r1)%A.s1); 
  A.set_v(-B.v%(A.s/A.s2));
}

/*
   Ideal inversion - Algorithm 0.17 of [B05]
   This function is a specialization of the ideal inversion function 
   for the case where s1 = 1, s2 = 1, [Ba] case. (i.e. H=1)
   
   A = B^(-1)

   -HZ
*/

void invV1a(cubic_ideal &A, cubic_ideal &B) {
  A.set_s(B.s);
  A.set_s1(B.s);   
  A.set_w(-B.u);
  A.set_v(-B.v);
}

/*
   Ideal inversion - Algorithm 0.18 of [B05]
   This function is a specialization of the ideal inversion function 
   for the case where s2 = 1. H != 1.

   A = B^(-1)
   
   -HZ
*/

void invV2(cubic_ideal &A, cubic_ideal &B) {
  ZZ_pX s_H = GCD(B.s,H);  
  ZZ_pX d, r1, r2; 

  A.set_s(B.s);
  A.set_s1(B.s/(B.s1*s_H));  
  A.set_s2(s_H);
 
  A.set_u((-B.w*H)%(B.s1*s_H)); 
  
  XGCD(d, r1, r2, H, A.s1); 
  
  A.set_w((-B.u*r1)%A.s1); 
  A.set_v(-((B.w*A.w*H)+B.v)%(A.s/A.s2)); 
}

/*
   Ideal inversion - Algorithm 0.19 of [B05]
   This function is a specialization of the ideal inversion function 
   for the case where s2 = 1, [Ba] case. (i.e. H=1)

   A = B^(-1)

   -HZ
*/

void invV2a(cubic_ideal &A, cubic_ideal &B) {
  A.set_s(B.s);
  A.set_s1(B.s/B.s1);   
  A.set_u(-B.w);
  A.set_w(-B.u);
  A.set_v(B.u*B.w-B.v);

  reduceCoeffs(A);
}

/*
   Ideal inversion - Algorithm 0.20 of [B05]
   Lemma 4.3.3 of [L09]
   This function performs ideal inversion for the general case. 
   
   A = B^(-1)

   -HZ
*/

void genInv(cubic_ideal &A, cubic_ideal &B) {
  ZZ_pX s_H = GCD(B.s,H); 
  ZZ_pX temp, d, r1, r2; 
  
  temp = B.s1*s_H;
  A.set_s(B.s);
  A.set_s1(B.s/temp);  
  A.set_s2(s_H/B.s2);
  
  A.set_u( -B.w*H%temp); 
  
  XGCD(d, r1, r2, H, A.s1); 
  A.set_w(-B.u*r1%A.s1);
  
  XGCD(d, r1, r2, B.s2, B.s/s_H); 
  A.set_v((-(B.v*B.s2*r1+B.w*A.w*H))%(A.s/A.s2)); 
}

/* 
   Ideal division - Algorithm 0.21 of [B05]
   This function is a specialization of the ideal division function 
   for the case where B.s1 = C.s1.

   -HZ
   A = C/B   - EJL
*/

void div1(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX Vtemp; 

  A.s = C.s/B.s;
  A.s1 = ZZ_pX(0,1);
  A.s2 = ZZ_pX(0,1);
  rem(A.u, C.u, A.s);

  Vtemp=C.v-C.w*A.u; 
  rem(A.v, Vtemp, A.s);
  A.w = ZZ_pX();

}


/* 
   Ideal division - Algorithm 0.21 of [B05]
   This function is a specialization of the ideal division function 
   for the case where s1 = 1 and s2 = 1.

   -HZ
   A = C/B   - EJL
*/

void divV1(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX temp;

  A.set_s(C.s/B.s);
  rem(temp, C.u, A.s);
  A.set_u(temp);
  rem(temp, C.v, A.s);
  A.set_v(temp);
  reduceCoeffs(A);
}

/* 
   Ideal division - Algorithm 0.22 of [B05]
   This function is a specialization of the ideal division function 
   for the case where there exist two primitive Ideals 
   
   -HZ

   A = C/B   - EJL

*/


void divV2a(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX H, U1, V1, temp, r1, r2, u, q, S_H;
  ZZ_pX d, d1, d2, d3;

  XGCD(d, r1, r2, (C.s/(C.s1*B.s2)), (B.s/(B.s1*B.s2)));
  GCD(d1, d, (B.u-C.u));
  GCD(d2, d1, (C.s2/B.s2));
  
  A.set_s((C.s*d2)/(B.s1*d1*B.s2));
  A.set_s1((C.s1*d1)/(B.s*d2));
  A.set_s2((C.s2)/(d2*B.s2));
  u = C.u+(C.w*H-B.u-C.u)*(r1*(C.s/(C.s1*d*B.s2)));

  XGCD(d3, r1, r2, 3*u*u, d/d1); 

  U1 = u - r1*((u*u*u+G)/d2); 
  rem(temp, U1, (A.s/A.s1));
  A.set_u(temp);

  rem(temp, C.w, A.s1);
  A.set_w(temp);
  q = (C.w-temp)/A.s1;

  GCD(S_H, A.s, H);
  XGCD(d, r1, r2, (A.s/S_H), (S_H/A.s2));
  V1 = (C.v-q*A.s1*A.u)*(r2*S_H/A.s2);
  rem(temp, V1, A.s/A.s2);
  A.set_v(temp);

  reduceCoeffs(A);

}

/* 
   Ideal division - Algorithm 0.23 of [B05]
   This function is a specialization of the ideal division function 
   for the case where there exist two primitive ideals with s2 = 1.

   -HZ
   A = C/B   - EJL
*/

void divV2b(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX U1, V1, r1, r2, u, q, S_H, Usqr;
  ZZ_pX d, d1, d2, d3;

  XGCD(d, r1, r2, (C.s/C.s1), (B.s/B.s1));

  GCD(d1, d, (B.u-C.u));
  
  A.set_s(C.s/(B.s1*d1));
  A.set_s1((C.s1*d1)/B.s);

  u = C.u+(C.w-B.u-C.u)*(r1*(C.s/(C.s1*d)));

  Usqr = sqr(u);

  XGCD(d2, r1, r2, 3*Usqr, d/d1);

  U1 = u - r1*((u*Usqr+G)/d2); 

  rem(A.u, U1, (A.s/A.s1));

  rem(A.w, C.w, A.s1);
  q = (C.w - A.w)/A.s1;

  V1 = (C.v-q*A.s1*A.u);

  rem(A.v, V1, A.s);

  reduceCoeffs(A);

}

/* 
   Ideal division - Algorithm 0.24 of [B05]
   This function is a specialization of the ideal division function 
   for the case where there exist two non-primitive ideals.

   -HZ

   Fixed.
   A = C/B  and F = GH^2 - EJL
*/

void divV3(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C){
  ZZ_pX d, d1, d2, r1, r2, u, q;

  // Step 1
  XGCD(d, r1, r2, C.s/C.s1, B.s/B.s1);

  // Step 2
  GCD(d1, d, B.u-C.u);
  
  // Step 3
  A.set_s(C.s/(d1*B.s1));
  A.set_s1(C.s1*d1/B.s);
  A.set_s2(ZZ_pX(0,1)); 

  // Step 4
  u = C.u + (C.w - B.u - C.u)*(C.s*r1/(C.s1*d));

  // Step 5
  XGCD(d2, r1, r2, d/d1, 3*sqr(u));

  // Step 6 and 7
  A.set_u((u-r2*((u*sqr(u) + f)/d2))%(A.s/A.s1));

  // Step 8
  A.set_w(C.w%A.s1);
  q = (C.w-A.w)/A.s1;

  // Step 9
  A.set_v((C.v-q*A.s1*A.u)%A.s);

}
/* 
   Ideal division - Algorithm 0.24 of [B05]
   This function is a specialization of the ideal division function 
   for the case where there exist two non-primitive ideals.

   -HZ

   Fixed.
   A = ((D)C)/B - EJL
*/

void divV4(ZZ_pX &D, cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  ZZ_pX D1, D2, Dtemp, S,S1;
  cubic_ideal temp, temp2; 

  temp2 = cubic_ideal(B);

  GCD(D1, B.s1, D);
  GCD(D2, B.s/(B.s1), D/D1);

  Dtemp = D1*D2;
  S=B.s/Dtemp; 
  S1=B.s1/D1; 

  ideal(temp, S, S1, temp2.u, temp2.v, temp2.w);

  temp2 = cubic_ideal(C);

  A = temp2/temp;

  if (!IsOne(Dtemp))
    {
      cubic_ideal I5;
      ideal(I5, Dtemp, D1, B.u, B.v, B.w);
      inverse_nored(temp, I5);

      A = A*temp;
    }
}
/*******************calculator functions************************/

  
/*******************overloaded operators************************/


/* 
   Overloaded multiplication operator for ideals (adapted from existing 
   code from Dr. Bauer)  
   
   -HZ

   Generalized to H != 1 -- EJL
*/

cubic_ideal cubic_ideal::operator* (const cubic_ideal& I2) const{
  cubic_ideal A, I1, result;
  ZZ_pX DD;
  I1.set_s(this->s);
  I1.set_s1(this->s1);
  I1.set_s2(this->s2);
  I1.set_u(this->u);
  I1.set_v(this->v);
  I1.set_w(this->w);

  GCD(DD, s, I2.s);

  if(IsOne(DD)) { // same as above
    
    if(s1 == I2.s1 == s2 == I2.s2 == 1) { // s' = s" = 1
      multV1(result, I1, I2);
    }
    else if(s2 == I2.s2 == 1) { // s" = 1, s != 1
      multV2(result, I1, I2);
    }
    
    else { 
      genMult(result, I1, I2);
    }
  }
  
  else {  // Not co-prime... either we're squaring or not.
    if(I1==I2){
      square(result, I1);
    }

    else if(IsOne(H)){ 
      ZZ_pX D1 = GCD3(I1.s/I1.s1, I2.s1, I1.u+I2.w);
      ZZ_pX D2 = GCD3(I2.s/I2.s1, I1.s1, I2.u+I1.w);
      ZZ_pX D3 = GCD3(I1.s1/D2, I2.s1/D1, I1.v + I2.v + I1.w*I2.w);
      
      ZZ_pX DD = D1*D2*D3;
      if(IsOne(DD)){
	mult30(result, I1, I2);
      }
      else{
	mult32(DD, D1, D2, D3, result, I1, I2);
      }
    }
    else{ 
      ZZ_pX s1H = GCD(I1.s, H), s2H = GCD(I2.s, H);
      ZZ_pX DH = GCD(I1.s2, I2.s2);
      ZZ_pX D1 = GCD3(I1.s/(I1.s1*s1H), I2.s1, I1.u+I2.w*H);
      ZZ_pX D2 = GCD3(I2.s/(I2.s1*s2H), I1.s1, I2.u+I1.w*H);
      ZZ_pX D3 = GCD3(I1.s1/D2, I2.s1/D1, I1.v + I2.v + I1.w*I2.w*H);
      ZZ_pX D4 = GCD(s1H/I1.s2, I2.s2);
      ZZ_pX D5 = GCD(s2H/I2.s2, I1.s2);
      
      ZZ_pX DD = DH*D1*D2*D3*D4*D5;
            
      if(IsOne(DD)){
	mult29(result, I1, I2);
      }
      else{
	mult31(DD, DH, D1, D2, D3, D4, D5, result, I1, I2);
      }
    }
  }
  reduce(A, result);

  return A;
}

/* 
   Overloaded division operator for ideals (adapted from existing 
   code from Dr. Bauer)  

   -HZ
*/

cubic_ideal  cubic_ideal::operator/(const cubic_ideal  &I2) {
  cubic_ideal I, I1;
  I1.set_s(this->s);
  I1.set_s1(this->s1);
  I1.set_s2(this->s2);
  I1.set_u(this->u);
  I1.set_v(this->v);
  I1.set_w(this->w);

  if(s1 == I2.s1 == s2 == I2.s2 == 1) { // s' = s" = 1

    if(I2.s % s == 0) {
      ZZ_pX temp1, temp2;
      rem(temp1, u, s);
      rem(temp2, I2.u, s);

      if(temp1 == temp2) {
        divV1(I, I2, I1);
	return I;
      }
    }
  }
  
  if(IsOne(s2) && IsOne(I2.s2)) { // s" = 1
    if(s1 == I2.s1){
      div1(I, I2, I1);
      return I;	
    }
    else {
      divV2b (I, I2, I1);
      return I;   
    }
  }
  divV2a(I, I2, I1);
  return I; 

}

/*
  Equals
*/

bool cubic_ideal::operator==(const cubic_ideal  & op2) const {
  cubic_ideal I;
  I.set_s(this->s);
  I.set_s1(this->s1);
  I.set_s2(this->s2);
  I.set_u(this->u);
  I.set_v(this->v);
  I.set_w(this->w);

  if((I.s == op2.s) && (I.s1 == op2.s1) &&(I.s2 == op2.s2) && (I.u == op2.u) &&(I.v == op2.v) && (I.w == op2.w))
    return true;
  else
    return false;
}

/* 
   Does not equal
*/

bool cubic_ideal::operator!=(const cubic_ideal  & op2) {
  cubic_ideal I;
  I.set_s(this->s);
  I.set_s1(this->s1);
  I.set_s2(this->s2);
  I.set_u(this->u);
  I.set_v(this->v);
  I.set_w(this->w);

  if(I.s != op2.s)
    return true;
  else if(I.s1 != op2.s1)
    return true;
  else if(I.s2 != op2.s2)
    return true;
  else if(I.u != op2.u)
    return true;
  else if(I.v != op2.v)
    return true;
  else if(I.w != op2.w)
    return true;
  else
    return false;
}


// Non-adjacent Form of n - [Reitswiesner, 1960]
// Input: an integer n
// Output: z = NAF(n), len = NumBits(z) 
// Algorithm 5.3.25 of [L09].
// -- EJL

void naf(ZZ n, int z[], int &len){
  len=0;
  while(n>0){
    if(n%2)
      z[len] = 2-to_int(n%4);
    else
      z[len] = 0;
    n-=z[len];
    n>>=1; //  n/=2;
    len++;
  }
}

/* 
   Overloaded operator for ideal exponentiation.
   Compute A = B^pow = B^k
   Uses the non-adjacent form (NAF) of pow, which is faster than
   standard binary exponentiation.
   -- EJL
*/

cubic_ideal  cubic_ideal::operator^(const ZZ& pow) const  {
  cubic_ideal A, B, C, X;
  ZZ k = pow;
  ZZ curPow = to_ZZ(2);
  int i, len;

  B.set_s(this->s);
  B.set_s1(this->s1);
  B.set_s2(this->s2);
  B.set_u(this->u);
  B.set_v(this->v);
  B.set_w(this->w);

  if(k==0)
    return C;
  else if(IsOne(k))
    return B;
  else if(k < 0){
    inverse(C,B);
    k = -k;
    return C^k;
  }
  else;

  A = cubic_ideal(B);
  inverse(X, B);
  
  // To store NAF(s). 
  int z[NumBits(k)+2];
  naf(k, z, len);

  for(i=len-2; i>=0; i--){
    C = A;
    square(A,C); // C = A^2
    if(z[i]==1){
      C = A;
      A = C*B;
    }
    else if(z[i]==-1){
      C = A;
      A = C*X;
    }
    else;
  }
  return A;
}

/*
   Ideal multiplication - A = B * C 
   This function simply forwards to the overloaded '*' function. 
   
   -HZ
*/

void cubic_ideal::multiply(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C) {
  A = B * C;
}

/* 
   Same as the * operator, except we don't reduce at the end.
   This is used for the (unit rank 0) reduction function since 
   the output is guaranteed to be reduced.
   Refer to Alg. 4.4.2 of [L09]

   A = B*C

*/

void multiply_nored(cubic_ideal  &A, const cubic_ideal  &B, const cubic_ideal  &C){
  ZZ_pX DD;

  GCD(DD, B.s, C.s);

  if(IsOne(DD)) { // same as above
    
    if(B.s1 == C.s1 == B.s2 == C.s2 == 1) { // s' = s" = 1
      multV1(A, B, C);
    }
    else if(B.s2 == C.s2 == 1) { // s" = 1, s != 1
      multV2(A, B, C);
    }
    
    else { 
      genMult(A, B, C);
    }
  }
  
  else {  // Not co-prime... either we're squaring or not.
    if(B==C){
      square(A, B);
    }

    else if(IsOne(H)){ 
      ZZ_pX D1 = GCD3(B.s/B.s1, C.s1, B.u+C.w);
      ZZ_pX D2 = GCD3(C.s/C.s1, B.s1, C.u+B.w);
      ZZ_pX D3 = GCD3(B.s1/D2, C.s1/D1, B.v + C.v + B.w*C.w);
      
      ZZ_pX DD = D1*D2*D3;
      if(IsOne(DD)){
	mult30(A, B, C);
      }
      else{
	mult32(DD, D1, D2, D3, A, B, C);
      }
    }
    else{ 
      ZZ_pX s1H = GCD(B.s, H), s2H = GCD(C.s, H);
      ZZ_pX DH = GCD(B.s2, C.s2);
      ZZ_pX D1 = GCD3(B.s/(B.s1*s1H), C.s1, B.u+C.w*H);
      ZZ_pX D2 = GCD3(C.s/(C.s1*s2H), B.s1, C.u+B.w*H);
      ZZ_pX D3 = GCD3(B.s1/D2, C.s1/D1, B.v + C.v + B.w*C.w*H);
      ZZ_pX D4 = GCD(s1H/B.s2, C.s2);
      ZZ_pX D5 = GCD(s2H/C.s2, B.s2);
      
      ZZ_pX DD = DH*D1*D2*D3*D4*D5;
            
      if(IsOne(DD)){
	mult29(A, B, C);
      }
      else{
	mult31(DD, DH, D1, D2, D3, D4, D5, A, B, C);
      }
    }
  }
}


// <D>A = B*C, Multiplication for infrastructure - Giant Steps
// Refer to Alg. 4.4.2 of [L09]
//  -- EJL

void multiply(cubic_ideal  &A, cubic_ideal  &B, cubic_ideal  &C, long &deg_diff){ 
  ZZ_pX D;

  GCD(D, B.s, C.s);

  if(IsOne(D)) { // same as above
    deg_diff = 0;
    if(B.s1 == C.s1 == B.s2 == C.s2 == 1) { // s' = s" = 1
      multV1(A, B, C);
    }
    else if(B.s2 == C.s2 == 1) { // s" = 1, s != 1
      multV2(A, B, C);
    }
    else { 
      genMult(A, B, C);
    }
  }
  
  else {  // Not co-prime... either we're squaring or not.
    if(B==C){
      square(A, B, deg_diff);
    }

    else if(IsOne(H)){ 
      ZZ_pX D1 = GCD3(B.s/B.s1, C.s1, B.u+C.w);
      ZZ_pX D2 = GCD3(C.s/C.s1, B.s1, C.u+B.w);
      ZZ_pX D3 = GCD3(B.s1/D2, C.s1/D1, B.v + C.v + B.w*C.w);
      
      D = D1*D2*D3;
      // Primitive.
      if(IsOne(D)){
	mult30(A, B, C);
	deg_diff = deg(A.s) - deg(B.s) - deg(C.s);
      }
      else{
	mult32(D, D1, D2, D3, A, B, C);
	deg_diff = deg(A.s) + deg(D) - deg(B.s) - deg(C.s);
      }
    }
    else{
      ZZ_pX s1H = GCD(B.s, H), s2H = GCD(C.s, H);
      ZZ_pX DH = GCD(B.s2, C.s2);
      ZZ_pX D1 = GCD3(B.s/(B.s1*s1H), C.s1, B.u+C.w*H);
      ZZ_pX D2 = GCD3(C.s/(C.s1*s2H), B.s1, C.u+B.w*H);
      ZZ_pX D3 = GCD3(B.s1/D2, C.s1/D1, B.v + C.v + B.w*C.w*H);
      ZZ_pX D4 = GCD(s1H/B.s2, C.s2);
      ZZ_pX D5 = GCD(s2H/C.s2, B.s2);
      
      D = DH*D1*D2*D3*D4*D5;
      
      if(IsOne(D)){
	mult29(A, B, C);
	deg_diff = deg(A.s) - deg(B.s) - deg(C.s);
      }
      else{
	mult31(D, DH, D1, D2, D3, D4, D5, A, B, C);
	deg_diff = deg(A.s) + deg(D) - deg(B.s) - deg(C.s);
      }
    }
  }
}

/*
   Ideal division - A = B / C : Only works for H=1.
   The algorithm says that division is performed for integral ideals 
   provided the result is also and integral ideal. Hence, we need a way 
   to check if the ideals and their result is also integral. 
   
   -HZ 
*/

void cubic_ideal::divide(cubic_ideal  & A,  cubic_ideal  & B, const cubic_ideal  & C) {
  A = B / C;
}

// Returns true if A==B and false otherwise.

bool cubic_ideal::IsEqual(const cubic_ideal &A, const cubic_ideal &B){
  return (A==B);
}

/*
   Ideal exponentiation - A = B ^ C 
  
*/

void cubic_ideal::power(cubic_ideal  & A, const cubic_ideal  & B, const ZZ & C) {
  A = B ^ C;
}
/*******************overloaded operators************************/


// Return 1 if order is probably the divisor class number of F_q(C).
// Return 0 if it is not.

int check(ZZ &order){
  cubic_ideal g, temp;

  for(int i = 0; i < 4; i++){
    random(g);
    temp = g^order;
    // This isn't the order. 
    if(!IsOne(temp.s)){
      return 0;
    }
  }
  return 1;
}

#endif /* CUBIC_IDEAL_CPP */
