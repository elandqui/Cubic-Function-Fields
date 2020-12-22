/*

  Contains utilities for computing an estimate, E, of the divisor class number,
  h, of the cubic function field K = Fq(C), C: Y^3 = f = GH^2. 

  Here, char(K) > 3.

 */
/* 
   hbounds.cc

   This file contains utilities for computing an estimate, E, of 
   the divisor class number, h, of the cubic function field K = Fq(C), 
   C: Y^3 = f = GH^2, and an upper bound, U, on the error |h-E|. 
   For the functions, we will use results and algorithms from:

   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)
   [SS7] R. Scheidler and A. Stein, "Class number approximation in 
                      cubic function fields"
   [SS8] R. Scheidler and A. Stein, "Approximating Euler products 
                      and class number computation in algebraic 
                      function fields"
   
   How to create the library:
   --------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download hbounds.cc, hbounds.h and makefile 
      into ../Hbounds.
   Compile libhbounds via make in ../Hbounds.

   Make sure that the paths to the NTL libraries are
   specified correctly in the makefile.
   
   Known Problems: 
   ---------------
   none
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   3.a. Possibly apply sieving in Sv1 and Sv2 for lambda > 2.
   3.b. rec3() rec32(), and isCube() can probably be improved via reciprocity.

   Author: 
   Eric Landquist
*/
#include "hbounds.h"

NTL_CLIENT

int factordegs[21]; // The degrees of the prime factors of G and H.

int lambda=0;    // The bound used for computing E and L.
double alpha=1.0;  // An optimizer to adjust L2;

long splitting[4];  // Stores the number of primes in each splitting category.
                    // 0: (P) = (p)^3 i.e. ramification.
                    // 1: (P) = (p)   i.e. inert
                    // 2: (P) = (p1)(p2)(p3) i.e. split completely
                    // 3: (P) = (p1)(p2) i.e. partial splitting
                    // (4:) (P) = (p1)(p2)^2 i.e. wild ramification
                    // Note: If q = 1 (mod 3), there is no partial splitting
                    // Note: If char(K) != 3, there is no wild ramification

vec_ZZ SvCache;  // Stores values of Sv(a) for 1 <= v <= lambda.
                 // If q = 1 mod 3, it stores Sv(1) and Sv(3).
                 // If q = 2 mod 3, it stores Sv(1) and Sv(2).


int MU[15] = {1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1};

int numPrimes = 320;
long primes[320] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823,  827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129};

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
ZZ U;

// Input: verbose - 0 for no output
//                - 1 for output
// Finds an approximation E of h and a bound U such that
// |h-E| < U.
void approxh(int verbose){
  int s1, s2;
  int n, v; // v|n
  int a, i, j, l=1;
  RR B;  // Set logE2 = A(K) + B.
  ZZ B1; // Used to find the trivial bound E1.
  RR C, S, T; // Temp variables.
  RR Q = to_RR(q);
  RR A;
  RR psi1, psi2, psi3;
  ZZ SvSum; // Used for psi3.
  ZZ zero = ZZ(), qm1 = q-1;

  /* Set lambda. See Eqn. 6.12 of [L09]. */

  if(genus%5 == 0)
    lambda += 2*genus/5;
  else if((genus%5 == 1)||(genus%5 == 4))
    lambda += (2*genus-2)/5;
  else if(genus%5 == 2)
    lambda += (2*genus-4)/5;
  else // if(genus%5 == 3)
    lambda += (2*genus-1)/5;

  if(lambda<=0)
    lambda=1;

  /**************************************/
  /* In here we store values of S_v(a). */
  /**************************************/

  // In the case q%3 == 2, the vector will be structured:
  // [S1(1) S1(2) S2(1) S2(2) S3(1) S3(2) ... ]
  // cacheFilled contains either a 0 or 1 depending on if Sv(a) has 
  // been computed or not.

  SvCache.SetLength(2*lambda);
  int cacheFilled[2*lambda]; 
  for(n = 0; n<2*lambda; n++)
    cacheFilled[n] = 0;

  if(verbose)
    cout<<"Running the approximation of h with LAMBDA = "<<lambda<<".\n\n";

  // First compute logE2 = A(K) + B + C. 
  // Begin with A(K) - The contribution of the infinite place(s).
  // See Eqn. 6.6 of [L09].
  if(deg(f)%3){
    s1=s2=0;
  } else if(q%3 == 2){
    s1=0; s2=-1;
  } else if(isCube(LeadCoeff(f))){
    s1=-2; s2=1;
  } else{
    s1=s2=1;
  }
  A = (genus+2)*log(Q) - log(sqr(Q)+to_RR(s1)*Q + to_RR(s2));
  logE1 = A;

  /***********************************************/
  /* Precompute the S_v(j) for 1 <= v <= lambda. */
  /***********************************************/

  // This will most likely be the longest running loop 
  // before the baby-step, giant-step or kangaroo portion(s).
  // See (6.6) of [L09] for the definition of S_v(j).
  Epr = exp(A);
  for(v=1; v<=lambda; v++){
    splitting[0] = splitting[1] = splitting[2] = splitting[3] = 0;
    if((q%3 == 1)||(v%2 == 0)) {
      Sv1(v, 1, zero, qm1);
      Sv1(v, 3, zero, qm1);
    } 
    else{
      Sv2(v, 1);
      Sv2(v, 2);
    }
    Epr*=(power(power(Q,2*v)/(power(Q,2*v) + power(Q,v) + 1.0), splitting[1])*power(power(Q,2*v)/(power(Q,2*v) -2.0*power(Q,v) + 1.0), splitting[2])*power(power(Q,2*v)/(power(Q,2*v) - 1.0), splitting[3]));
  }
  Ep = RoundToZZ(Epr);
  /*************************/
  /* Approximate h with E. */
  /*************************/  
  // See Sections 6.3.4 and 6.3.6 of [L09] for the derivation and calulation
  // of the three estimates.
  // Now we compute C = SUM_{n=1}^{oo} 1/nq^n SUM_{v|n, v <= lambda} vS_v(n/v).
  C = RR();

  B = to_RR(1);
  for(i=1; i <= lambda; i++){
    B *= Q; //power(Q,i);
        
    // This is just to check the trivial bounds and 
    // is included to empirically show that the better bounds 
    // are better. We only include this now for the case q=1%3.
    B1 = ZZ();
    for(j=1; j<=(lambda/i); j++){
      if((q%3 == 1) || ((i%2) == 0)){
	if(j%3 == 0){
	  B1 +=  SvCache[2*i-1];
	}	  
	else{
	  B1 += SvCache[2*i-2];
	}	
      }
      else{
	if(j%2 == 0){
	  B1 +=  SvCache[2*i-1];
	}	  
	else{
	  B1 += SvCache[2*i-2];
	}	
      }
      logE1 += to_RR(B1)/to_RR(j*power(B, j));
    }
    // This is to compute E2.
    if((q%3 == 1) || ((i%2) == 0)){
      S = power(B,3);
      C += (to_RR(SvCache[2*(i-1)] - SvCache[2*i-1])*log1p(-inv(S))/3.0 - to_RR(SvCache[2*(i-1)])*log1p(-inv(B)));
    }
    else{
      S = sqr(B);
      C+= (to_RR(SvCache[2*(i-1)] - SvCache[2*i-1])*log1p(-inv(S))/2.0 - to_RR(SvCache[2*(i-1)])*log1p(-inv(B)));
    }
  }
  logE2 = A + C;
  E2 = exp(logE2);
  //E2 = Ep;
  // This is the trivial bound.
  E1 = RoundToZZ(exp(logE1));
  RoundToZZ(E, E2);

  /************************************************/
  /* Approximate the square root of the error, L. */
  /************************************************/

  // First, B = SUM_{n=lambda+1}^{oo} 1/nq^n SUM_{v|n, v > lambda} vS_v(n/v).
  psi2=RR();

  // The first step is estimating the first lambda terms of this:
  // S_{lambda+1}(1)/(q^{lambda+1}) + ... + S_{2*lambda}(1)/(q^{2*lambda})
  RR Tinv = inv(sqrt(Q));
  T = sqrt(Q);
  
  /*********************************************/
  /* These are the trivial bounds: E1 and L1^2 */
  /*********************************************/
  
  E1 = RoundToZZ(exp(logE1));
  
  psi1 = -log1p(-Tinv); 
  for(i=1; i<=lambda; i++){
    psi1 -= power(Tinv,i)/i;
  }
  psi1*=(2*genus);
  psi1 -= 2*(log1p(-inv(Q)));
  for(i=1; i<=lambda; i++){
    psi1 -= (2*inv(i*power(Q,i)));
  }
  
  L12 = CeilToZZ(exp(logE1)*expm1(psi1)+1.0/2.0);
  
  /*************************************/
  /* These are the bounds: E2 and L2^2 */
  /*************************************/
  
  i = lambda+1;
  n = i+1;
  for(a=0; a<15; a++){
    if(i%primes[a] == 0){
      l = primes[a];
      break;
    }
  }
  psi2 = (2.0/(i*power(Q,i)))*(1.0+(Q/(Q-1.0))*(power(Q,i/l)-1.0)) + 2*genus*power(Tinv, i)/i + 2*genus*T*power(Tinv, n)/(n*(T-1.0)) + 4*Q*pow(Q, (to_RR(l)-1.0)/to_RR(l))*pow(Q,-to_RR(n*(l-1))/to_RR(l))/(n*(Q-1.0)*(pow(Q, (to_RR(l)-1.0)/to_RR(l))-1.0));

  // We're taking advantage of easy inversions here.
  L2 = CeilToZZ(E2*expm1(psi2)+1.0/2.0);
  L = CeilToZZ(sqrt(alpha*to_RR(L2)));
  
  // In the off-chance that an estimate is partly outside
  // of the Hasse interval, we correct it.
  // (sqrt(q) - 1)^(2genus) < h < (sqrt(q) + 1)^(2genus).
  A = power(T - 1, 2*genus) - to_RR(E - L2);
  if(A>0){
    L2 = E - CeilToZZ(power(T - 1, 2*genus));
    L = CeilToZZ(sqrt(alpha*to_RR(L2)));
  }
  A = to_RR(E + L2) - power(T + 1,2*genus) ;
  if(A>0){
    L2 = FloorToZZ(power(T + 1, 2*genus)) - E;
    L = CeilToZZ(sqrt(alpha*to_RR(L2)));
  }
  
    
  /*************************************/
  /* These are the bounds: E3 and L3^2 */
  /*************************************/
  
  // E3 = E2.
  
  // This is a little bit of a speed-up.
  // Instead of dividing odd i by 3,
  // and even i by 2,
  // we divide i by its smallest prime factor.
  // There isn't a speed-up until i=5 is reached,
  // i.e. unless lambda > 2.
  
  for(a=0; a<15; a++){
    if(i%primes[a] == 0){
      l = primes[a];
      break;
    }
  }
  S = 1.0/to_RR(l);
  C = pow(Q, S);
  S = Q/C;
  
  // Finding |SUM(x_j^{lambda+1})|
  if(deg(f)%3)
    psi3 = 2.0;
  else
    psi3 = RR();
  
  SvSum = ZZ();
  // Finding |SUM_{v|(lambda+1), v != lambda +1}vS_v((lambda+1)/v)|
  for(j=1; j<i; j++){
    if(i%j == 0){
      if((q%3 == 1) || ((j%2) == 0)){
	if((i/j)%3)
	  SvSum += j*SvCache[2*(j-1)];
	else
	  SvSum += j*SvCache[2*j-1];
      }
      else{
	if((i/j)%2)
	  SvSum += j*SvCache[2*(j-1)];
	else
	  SvSum += j*SvCache[2*j-1];
      }
    }
  }

  psi3 += to_RR(abs(SvSum));
  psi3 /= (i*power(Q, i));
  
  psi3 +=  2*genus*power(Tinv, i)/i + 2*genus*T*power(Tinv, n)/(n*(T-1.0)) + 4*Q*pow(Q, (to_RR(l)-1.0)/to_RR(l))*pow(Q,-to_RR(n*(l-1))/to_RR(l))/(n*(Q-1.0)*(pow(Q, (to_RR(l)-1.0)/to_RR(l))-1.0));

  U = L32 = CeilToZZ(E2*expm1(psi3)+1.0/2.0);
  
  // I'm going with this guy to speed things up.
  L = CeilToZZ(sqrt(alpha*to_RR(L32)));
  
  // Previously, there was a correction if the new interval
  // went out of the bounds of the Hasse interval.
  // That approach changed E, so we won't do that anymore.
  
  if(verbose){

    cout<<"Hasse Interval: ["<<CeilToZZ(power(T - 1, 2*genus))<<", "<<FloorToZZ(power(T + 1, 2*genus))<<"].\n";
    //cout<<"New interval 1: ["<<E1-L12<<", "<<E1+L12<<"].\n";
    //cout<<"New interval 2: ["<<E-L2<<", "<<E+L2<<"].\n";
    cout<<"New interval 3: ["<<E-U<<", "<<E+U<<"].\n\n";
    //cout<<"Estimate of the class number E_1 = "<<E1<<".\n";
    cout<<"Estimate of h:                 E = "<<E<<".\n";
    //cout<<"Our bounds on the error:   L_1^2 = "<<L12<<".\n";
    //cout<<"Our bounds on the error:   L_2^2 = "<<L2<<".\n";
    cout<<"Upper bound on the error:      U = "<<U<<".\n\n";
    //cout<<"Setting L = "<<L<<".\n\n";
  }
}

// Computes the value S_v(a) = SUM_{deg(p)=v}(z1(p)^a + z2(p)^a)
// in the case that q^v = 1 (mod 3).
// See Alg. 6.3.10 of [L09] see also Sections 6.3.4 and 6.3.6 of [L09].
// Input: v - nu in [L09]
//        a - input of S_v(a).
//        start - beginning of polynomial block if parallelized.
//        end - end of polynomial block if parallelized.
// Output: stores function value in SvCache.
// Array structure: 
// q = 1 (mod 3): [S_1(1), S_1(3), S_2(1), S_2(3), ...]
// q = 2 (mod 3): [S_1(1), S_1(2), S_2(1), S_2(3), ...]
void Sv1(int v, int a, ZZ &start, ZZ &end){
  int loc = 2*(v-1) + (a%3==0), i, divlim;
  ZZ T=ZZ(), j, c;
  ZZ_pX P;

  // In this case, z1(P)^a + z2(P)^a = 2 if P!|GH and 0 if P|GH.
  if(a%3 == 0){
    divlim = v/4;
    T = ZZ();
    for(i=1; i<=divlim; i++){
      if(v%i == 0)
	T+=mu(v/i)*power(q,i);
    }
    if(v%3==0)
      T-=power(q,v/3);
    if(v%2==0)
      T-=power(q,v/2);
    T+=power(q,v);

    T/=v;

    // Subtract off the number of irreducible factors of GH of degree v.
    T-=factordegs[v-1];
    
    SvCache[loc] = 2*T;
    return;
  }

  // Otherwise,  z1(P)^a + z2(P)^a = 2 if (f|P)_3 = 1 and 
  // z1(P)^a + z2(P)^a = -1 if (f|P)_3 != 1
  P = ZZ_pX(v, 1);
  j = start; //j=ZZ();
  splitting[3] = 0;
  // Loop through every monic irreducible degree v polynomial.
  // In the case that deg(P) = v = 1, every polynomial is irred.
  if(v == 1){
    splitting[0] = factordegs[0];
    // If there are some ramified primes of this degree, 
    // We need to check those.
    if(factordegs[0]){
      while(j <= end){
	// Get the next polynomial.
	SetCoeff(P, 0, to_ZZ_p(j));
	if(IsZero(f%P)); // This test can be sped up.
	else if(chi(P) == 1){
	  // This is the complete splitting case.
	  splitting[2]++;
	  T+=2;
	}
	else{
	  // This is the inert case.
	  splitting[1]++;
	  T--;
	}
	j++;
      }
    }
    else{
      while(j <= end){
	// Get the next polynomial.
	SetCoeff(P, 0, to_ZZ_p(j));
	if(chi(P) == 1){
	  // This is the complete splitting case.
	  splitting[2]++;
	  T+=2;
	}
	else{
	 // This is the inert case.
	  splitting[1]++;
	  T--;
	}
	j++;
      }
    }
  }
  // We use a nice trick to run through every irreducible 
  // degree 2 polynomial.
  else if (v == 2){
    splitting[0] = factordegs[1];
    if(IsZero(start)) start++;
    //P = ZZ_pX(v, 1);
    // If there are some ramified primes of this degree, 
    // We need to check those.
    if(factordegs[1]){
      if((q%3) == 1){
	for(j = to_ZZ(1); j <= end; j++){
	  // If x^2 - j is irreducible, cycle through each (x+c)^2 - j.
	  if(Jacobi(j, q) == -1){
	    SetCoeff(P, 0, to_ZZ_p(-j));
	    for(c = ZZ(); c < q; c++){
	      if(IsZero(f%P)); // This test can be sped up.
	      else if(chi(P) == 1){
		// This is the complete splitting case.
		splitting[2]++;
		T+=2;
	      }
	      else{
		// This is the inert case.
		splitting[1]++;
		T--;
	      }
	      // Get the next polynomial.
	      SetCoeff(P, 0, to_ZZ_p(rep(coeff(P,0))+(2*c+1)));
	      SetCoeff(P, 1, to_ZZ_p(rep(coeff(P,1))+2));
	    }
	  }
	}
      }
      else{ // q = 2 (mod 3)
	ZZ Q = (sqr(q)-1)/3;
	for(j = to_ZZ(1); j <= end; j++){
	  // If x^2 - j is irreducible, cycle through each (x+c)^2 - j.
	  if(Jacobi(j, q) == -1){
	    SetCoeff(P, 0, to_ZZ_p(-j));
	    //cout<<P<<endl;
	    for(c = ZZ(); c < q; c++){
	      if(IsZero(f%P)); // This test can be sped up.
	      //cout<<P<<" "<<chi2(P, Q)<<endl;
	      else if(IsOne(chi2(P, Q))){
		// This is the complete splitting case.
		splitting[2]++;
		T+=2;
	      }
	      else{
		// This is the inert case.
		splitting[1]++;
		T--;
	      }
	      // Get the next polynomial.
	      SetCoeff(P, 0, to_ZZ_p(rep(coeff(P,0))+(2*c+1)));
	      SetCoeff(P, 1, to_ZZ_p(rep(coeff(P,1))+2));
	    }
	  }
	}
      }
    }
    else{
      if((q%3) == 1){
	for(j = to_ZZ(1); j <= end; j++){
	  // If x^2 - j is irreducible, cycle through each (x+c)^2 - j.
	  if(Jacobi(j, q) == -1){
	    SetCoeff(P, 0, to_ZZ_p(-j));
	    //cout<<P<<endl;
	    for(c = ZZ(); c < q; c++){
	      if(chi(P) == 1){
		// This is the complete splitting case.
		splitting[2]++;
		T+=2;
	      }
	      else{
		// This is the inert case.
		splitting[1]++;
		T--;
	      }
	      // Get the next polynomial.
	      SetCoeff(P, 0, to_ZZ_p(rep(coeff(P,0))+(2*c+1)));
	      SetCoeff(P, 1, to_ZZ_p(rep(coeff(P,1))+2));
	    }
	  }
	}
      }
      else{ // q = 2 (mod 3)
	ZZ Q = (sqr(q)-1)/3;
	for(j = to_ZZ(1); j <= end; j++){
	  // If x^2 - j is irreducible, cycle through each (x+c)^2 - j.
	  if(Jacobi(j, q) == -1){
	    SetCoeff(P, 0, to_ZZ_p(-j));
	    //cout<<P<<endl;
	    for(c = ZZ(); c < q; c++){
	      //cout<<P<<" "<<chi2(P, Q)<<endl;
	      if(IsOne(chi2(P, Q))){
		// This is the complete splitting case.
		splitting[2]++;
		T+=2;
	      }
	      else{
		// This is the inert case.
		splitting[1]++;
		T--;
	      }
	      // Get the next polynomial.
	      SetCoeff(P, 0, to_ZZ_p(rep(coeff(P,0))+(2*c+1)));
	      SetCoeff(P, 1, to_ZZ_p(rep(coeff(P,1))+2));
	    }
	  }
	}
      }
    }
  }
  else{
    SetCoeff(P, v-1, to_ZZ_p(start));
    // Here we have deg(P) > 2, so we test for irreducibility.
    splitting[0] = factordegs[v-1];
    if(factordegs[v-1]){
      if((q%3) == 1){
	while(IsOne(LeadCoeff(P)) && (rep(coeff(P, lambda-1)) <= end)){
	  // Formerly if(ProbIrredTest(P, iter=1)){
	  if(DetIrredTest(P)){
	    if(IsZero(f%P));
	    else if(chi(P) == 1){
	      splitting[2]++;
	      T+=2;
	    }
	    else{
	      splitting[1]++;
	      T--;
	    }
	  }
	}
	// Get the next polynomial.
	SetCoeff(P, 0, coeff(P,0)+1);
	i=0;
	// "Carry the 1."
	while(coeff(P,i)==0){
	  i++;
	  SetCoeff(P, i, coeff(P,i)+1);
	}
      }
      else{ // q = 2 (mod 3)
	ZZ Q = (power(q, v)-1)/3;
	while(IsOne(LeadCoeff(P)) && (rep(coeff(P, lambda-1)) <= end)){
	  // Formerly if(ProbIrredTest(P, iter=1)){
	  if(DetIrredTest(P)){
	    if(IsZero(f%P));
	    else if(chi2(P, Q) == 1){
	      splitting[2]++;
	      T+=2;
	    }
	    else{
	      splitting[1]++;
	      T--;
	    }
	  }
	}
	// Get the next polynomial.
	SetCoeff(P, 0, coeff(P,0)+1);
	i=0;
	// "Carry the 1."
	while(coeff(P,i)==0){
	  i++;
	  SetCoeff(P, i, coeff(P,i)+1);
	}
      }
    }
    else{
      if((q%3) == 1){
	while(IsOne(LeadCoeff(P)) && (rep(coeff(P, lambda-1)) <= end)){
	  
	  if(DetIrredTest(P)){
	    if(chi(P) == 1){
	      splitting[2]++;
	      T+=2;
	    }
	    else{
	      splitting[1]++;
	      T--;
	    }
	  }
	  // Get the next polynomial.
	  SetCoeff(P, 0, coeff(P,0)+1);
	  i=0;
	  // "Carry the 1."
	  while(coeff(P,i)==0){
	    i++;
	    SetCoeff(P, i, coeff(P,i)+1);
	  }
	}
      }
      else{ // q = 2 (mod 3)
	ZZ Q = (power(q, v)-1)/3;
	while(IsOne(LeadCoeff(P)) && (rep(coeff(P, lambda-1)) <= end)){
	  // Formerly if(ProbIrredTest(P, iter=1)){
	  if(DetIrredTest(P)){
	    if(IsZero(f%P));
	    else if(chi2(P, Q) == 1){
	      splitting[2]++;
	      T+=2;
	    }
	    else{
	      splitting[1]++;
	      T--;
	    }
	  }
	}
	// Get the next polynomial.
	SetCoeff(P, 0, coeff(P,0)+1);
	i=0;
	// "Carry the 1."
	while(coeff(P,i)==0){
	  i++;
	  SetCoeff(P, i, coeff(P,i)+1);
	}
      }
    }
  }
  
  SvCache[loc] = T;
  return;

}

// Computes the value S_v(a) = SUM_{deg(p)=v}(z1(p)^a + z2(p)^a)
// in the case that q^v = 2 (mod 3)
// See Alg. 6.3.11 of [L09] see also Sections 6.3.4 and 6.3.6 of [L09].
// Input: v - nu in [L09]
//        a - input of S_v(a).
// Output: stores function value in SvCache.
// Array structure: 
// q = 2 (mod 3): [S_1(1), S_1(2), S_2(1), S_2(3), ...]
void Sv2(int v, int a){
  int loc = 2*(v-1) + (a%2==0), i, divlim;
  ZZ T;

  // If v is odd and a is odd, then z1(p)^a + z2(p)^a = 0.

  if((v%2) && (a%2)){
    SvCache[loc] = ZZ();
    return;
  }
 
  // If v is even, or v is odd and a is even, then z1(p)^a + z2(p)^a = 2.
  // So we just need to calculate the number Lq(n) of monic irreducible
  // polynomials of degree n over Fq, subtract the number of irreducible
  // factors of G and H of degree v, multiply by 2, and return.
  divlim = v/4;
  T = ZZ();
  for(i=1; i<=divlim; i++){
    if(v%i == 0)
      T+=mu(v/i)*power(q,i);
  }
  if(v%3==0)
    T-=power(q,v/3);
  if(v%2==0)
    T-=power(q,v/2);
  T+=power(q,v);

  T/=v;

  // Subtract off the number of irreducible factors of GH of degree v.
  T-=factordegs[v-1];
  splitting[0] = factordegs[v-1];
  splitting[1] = splitting[2] = 0;
  splitting[3] = to_long(T);

  T*=2;
  SvCache[loc] = T;
  return;
}

// The Moebius mu function.
// See any text in analytic number theory for its definition.
int mu(int a){
  int i, k=0;

  if(a<16)
    return(MU[a-1]);

  // This will only fail if there is a square factor that is not 47-smooth
  // Or if there are more than 2 large primes.

  for(i=0; (i<15) && (a!=1); i++){
    if(a%i == 0){
      // If a is not squarefree, return 0.
      if((a/i)%i == 0)
	return(0);
      else
	k++;
    }
  }
  if(a==1) {
    if(k%2) return(1);
    else return(-1);
  }  else if(ProbPrime(a)){
    if(k%2) return(-1);
    else return(1);
  }  else{
    // Here, I'm assuming that the composite remainder has
    // two prime factors. It will also work if there are
    // an even number of large primes.
    if(k%2) return(1);
    else return(-1);
  }

}

// Computes the character [f|P]_3 if q = 1 (mod 3).
// Follows Algorithm 6.2 of [SS7] and Alg. 6.3.7 of [L09].
ZZ_p chi(const ZZ_pX &_P){
  int df, p;
  ZZ_p e = to_ZZ_p(1);
  ZZ_pX F = f, T, P = _P;

  p = deg(P)%3;
  while(deg(F)>0){
    F%=P;
    df = deg(F)%3;
    e*=(power(rec3(LeadCoeff(F)), p)*power(inv(rec3(LeadCoeff(P))),df));
    T = F;
    F = P;
    P = T;
    p = df;
  }  
  return(e*power(rec3(LeadCoeff(F)), deg(P)%3));

}

// Computes the character [f|P]_3 if q = 2 (mod 3).
// See Eqn. 6.13 of [L09].
ZZ_pX chi2(const ZZ_pX &P, const ZZ &Q){
  const ZZ_pX F = f%P;
  ZZ_pXModulus T(P);
  return PowerMod(F, Q, T);
}

// Computes the cubic reciprocity of a mod q.
// rec3(): Input: a - an element of F_q.
//         Output: 1 if a is a cube in F_q
//                 something else if not.
// rec32(): Input: a - an element of F_Q -- an extension of F_q.
//                 Q - a power of q.
//          Output: 1 if a is a cube in F_Q
//                  something else if not.
inline ZZ_p rec3(ZZ_p a){
  return(power(a, (q-1)/3));
}
inline ZZ_p rec32(ZZ_p a, ZZ &Q){
  return(power(a, (Q-1)/3));
}

// Returns h = a (mod l)
// Only applies to rank 0 curves.

void smallOrder0(ZZ &a, ZZ &l){
  ZZ k=to_ZZ(3); 
  l = to_ZZ(1);
  vec_pair_ZZ_pX_long factorsG, factorsH;

  if((rank > 0) && (genus > 3)){
    a = ZZ();
    l = to_ZZ(1);
    return;
  }

  factorsG = berlekamp(inv(LeadCoeff(G))*G);
  
  if(IsOne(H)){
    power(l, k, factorsG.length()-1);
  }
  else{
    factorsH = berlekamp(inv(LeadCoeff(H))*H);
    power(l, k, factorsG.length()+factorsH.length()-2);
    if(IsOne(l)){
      a = ZZ();
      if(deg(H) < genus - 1)
	l = to_ZZ(3);
      return;
    }      
  }
  if(!IsOne(l))
    a = ZZ();
  else if(q%3 == 1){
    a = to_ZZ(1);
    l = to_ZZ(3);
  }
  else
    a = ZZ();
    
  return;
}

// Checks to see if ord = 1 mod 3 is a possible order
// returns 1 if it is a possible order
// returns 0 if it is not.
// This function really only applies to really small examples.
int checkPrimes(ZZ &ord){
  ZZ test = ord;
  int k=0;
  int i;

  if(rank == 0){
    if(deg(H) != 0)
      return 1;
  }
  else if(rank == 1)
    return 1;
  else{
    //if(!((deg(G) == 4) && (deg(H) == 1)))
    if(q%2 == 3)
      return 1;
  }

  while(test%3 == 0){
    test/=3;
  }

  if(test%3 == 2)
    return 0;

  // Check all the primes up to 2129.
  for(i=0; (i<numPrimes) && !IsOne(test); i++){
    k=0;
    while(test%primes[i] == 0){
      test/=primes[i];
      k++;
    }
    if(primes[i]%3 == 2){
      if(k%2)
	return 0;
    }
  }

  if(test%3 == 1){
    return 1;
  }

  return 0;

  // If it gets here, then factor has a couple factors larger than 2129,
  // an odd number of which are = 2 mod 3. We should improve this to pull
  // out more factors    
}

// Check an order in the unit rank 0 case.
// Input: order - a potential class number
//        factor - a known factor of the class number
// Output: 1 (true) if order is a likely order based on factor
//         0 (false) if order/factor contains primes not dividing factor.

int shareFactors(ZZ order, ZZ factor){
  ZZ mult = order/factor;
  ZZ div;

  if(IsZero(factor % mult))
    return 1;

  div = GCD(mult, factor);

  if(IsOne(div))
    return 0;

  if(IsZero(factor % (mult/div)))
    return 1;

  // Strip small primes common to factor and mult and compare again.
  for(int i=0; i<4; i++){
    if(div % primes[i]){
      while ( div%primes[i] == 0){
	div/=primes[i];
	mult/=primes[i];
      }
      while (mult%primes[i] == 0){
	mult/=primes[i];
      }
    }
  }

  if(IsZero(factor % (mult/div)))
    return 1;

  return 0;
  
}


// Returns 1 if a is a cube in ZZ_p.
// Returns 0 if not.

int isCube(ZZ_p a){
  if(q%3 == 2)
    return(1);
  else
    return(IsOne(power(a, (q-1)/3))); 
}
