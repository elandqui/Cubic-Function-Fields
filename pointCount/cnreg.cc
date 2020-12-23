#include "cnreg.h"

NTL_CLIENT

// Finds an approximation E of h and a bound L such that
// |h-E|<L^2.

void approxh(){
  int s1, s2;
  int n, v, divlim; // v|n
  int a, i, j, l;
  RR B;  // Set logE2 = A(K) + B.
  ZZ B1; // Used to find the trivial bound E1.
  RR C, S, T; // Temp variables.
  RR Q = to_RR(q);
  RR A;
  RR psi1, psi2, psi3;
  ZZ SvSum; // Used for psi3.

  //RR::SetPrecision(300);

  /***************/
  /* Set lambda. */
  /***************/

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

  /*************************************/
  /* In here we store values of Sv(a). */
  /*************************************/

  // In the case q%3 == 2, the vector will be structured:
  // [S1(1) S1(2) S2(1) S2(2) S3(1) S3(2) ... ]
  // cacheFilled contains either a 0 or 1 depending on if Sv(a) has 
  // been computed or not.

  SvCache.SetLength(2*lambda);
  int cacheFilled[2*lambda]; 
  for(n = 0; n<2*lambda; n++)
    cacheFilled[n] = 0;

  cout<<"Running the approximation of h with LAMBDA = "<<lambda<<".\n\n";

  // First compute logE2 = A(K) + B + C. 
  // Begin with A(K).
  if(deg(D)%3){
    s1=s2=0;
  } else if(q%3 == 2){
    s1=0; s2=-1;
  } else if(isCube(LeadCoeff(D))){
    s1=-2; s2=1;
  } else{
    s1=s2=1;
  }
  A = (genus+2)*log(Q) - log(sqr(Q)+to_RR(s1)*Q + to_RR(s2));
  logE1 = A;
  //cout<<"So far our guess for h is A(K) = "<<to_ZZ(exp(A))<<".\n";
 
  /***********************************************/
  /* Precompute the S_v(j) for 1 <= v <= lambda. */
  /***********************************************/

  // This will most likely be the longest running loop 
  // before the baby-step, giant-step portion.
  Epr = exp(A);
  for(v=1; v<=lambda; v++){
    splitting[0] = splitting[1] = splitting[2] = 0;
    if((q%3 == 1)||(v%2 == 0)) {
      Sv1(v, 1);
      Sv1(v, 3);
    } 
    else{
      // Note: This must be fixed for q=2(mod 3)
      // since S_v(a) takes on four values in this case:
      // S_v(1), S_v(2), S_v(3), and S_v(6).
      Sv2(v, 1);
      Sv2(v, 2);
    }
    //cout<<splitting[0]<<" "<<splitting[1]<<" "<<splitting[2]<<endl;
    Epr*=(power(power(Q,2*v)/(power(Q,2*v) + power(Q,v) + 1.0),splitting[1])*power(power(Q,2*v)/(power(Q,2*v) -2.0*power(Q,v) + 1.0),splitting[2]));
  }

  Ep = RoundToZZ(Epr);

  /*************************/
  /* Approximate h with E. */
  /*************************/  

  // Now we compute C = SUM_{n=1}^{oo} 1/nq^n SUM_{v|n, v <= lambda} vS_v(n/v).
  C = RR();
  if(q%3 == 1){
    B = to_RR(1);
    for(i=1; i <= lambda; i++){
      B *= Q; //power(Q,i);
      S = power(B,3);

      // This is just to check the trivial bounds and 
      // is included to empirically show that the better bounds 
      // are better. We only include this now for the more
      // interesting case q=1%3.
      B1 = ZZ();
      for(j=1; j<=i; j++){
	if((i%j) == 0){
	  // B1 += v*S_v(n/v)
	  if((i/j)%3 == 0){
	    // B1 += v*S_v(3)
	    B1 +=  j*SvCache[2*j-1];
	  }	  
	  else{
	    // B1 += v*S_v(1)
	    B1 += j*SvCache[2*j-2];
	  }	  
	}
      }
      logE1 += to_RR(B1)/to_RR(i*B);

      // This is to compute E2.
      // logE2 = A(K) + SUM_{i=1}^{lambda} 
      //         (1/3)(S_v(1)-S_v(3))ln(1-1/q^(3v)) - S_v(1)ln(1-1/q^v)
      //C += (to_RR(SvCache[2*(i-1)] - SvCache[2*i-1])*log(1.0-inv(S))/3.0 - to_RR(SvCache[2*(i-1)])*log(1.0-inv(B)));      
      C += (to_RR(SvCache[2*(i-1)] - SvCache[2*i-1])*log1p(-inv(S))/3.0 - to_RR(SvCache[2*(i-1)])*log1p(-inv(B)));
    }
  }
  else {
    for(i=1; i <= lambda; i++){
      if(i%2){	
	B = power(Q,2*i);
	C+= to_RR(SvCache[2*i-1])*log(1.0-inv(B))/2;	
      }
      else {
	B = power(Q,i);
	S = sqr(B);
	C+= (to_RR(SvCache[2*(i-1)])*log((B+1.0)/(B-1.0)) - to_RR(SvCache[2*i-1])*log(1.0-inv(S)))/2;
      }
    }
  }
  
  logE2 = A + C;
  E2 = exp(logE2);
  // This is the trivial bound.
  E1 = RoundToZZ(exp(logE1));
  //cout<<"So far our guess for h is   E2 = "<<E2<<".\n";
  RoundToZZ(E, E2);
  //cout<<"Setting                      E = "<<E<<".\n\n";
  //cout<<"Testing ................... E' = "<<Ep<<".\n";
  //cout<<"Estimate of h:               E = "<<E<<".\n\n";

  /************************************************/
  /* Approximate the square root of the error, L. */
  /************************************************/

  // First, B = SUM_{n=lambda+1}^{oo} 1/nq^n SUM_{v|n, v > lambda} vS_v(n/v).
  psi2=RR();

  // The first step is estimating the first lambda terms of this:
  // S_{lambda+1}(1)/(q^{lambda+1}) + ... + S_{2*lambda}(1)/(q^{2*lambda})
  //if(q%3 == 1){
    //C = power(Q,lambda);
    T = sqrt(Q);

    /*********************************************/
    /* These are the trivial bounds: E1 and L1^2 */
    /*********************************************/

    E1 = RoundToZZ(exp(logE1));

    psi1 = log(T) - log(T-1.0);
    for(i=1; i<=lambda; i++){
      psi1 -= inv(i*power(T,i));
    }
    psi1*=(2*genus);
    psi1 += 2*(log(Q) - log(Q-1));
    for(i=1; i<=lambda; i++){
      psi1 -= (2*inv(i*power(Q,i)));
    }

    L12 = CeilToZZ(exp(logE1)*expm1(psi1)+1.0/2.0);

    /*************************************/
    /* These are the bounds: E2 and L2^2 */
    /*************************************/

    i = lambda+1;
    n = i+1;

    psi2 = 2.0*genus/(i*power(T,i)) + ((2*genus + 4)*T*Q*power(T, -n))/(n*(T-1.0)*(Q-1.0)) + 2.0/(i*power(Q, i));

    if(i%2)
      psi2 += (2.0*(Q/(Q-1.0))*power(Q, -i)*(power(Q,i/3)-1.0)/i);
    else
      psi2 += (2.0*(Q/(Q-1.0))*power(Q, -i)*(power(Q,i/2)-1.0)/i);	
      
    //L = CeilToZZ(sqrt(2*E2*expm1(B)+1.0));
    // We're taking advantage of easy inversions here.
    L2 = CeilToZZ(E2*expm1(psi2)+1.0/2.0);
    L = CeilToZZ(sqrt(alpha*to_RR(L2)));

    // In the off-chance that an estimate is partly outside
    // of the Hasse interval, we correct it.
    // (sqrt(q) - 1)^(2genus) < h < (sqrt(q) + 1)^(2genus).
    A = power(T - 1, 2*genus) - to_RR(E - L2);
    if(A>0){
      L2 = E - CeilToZZ(power(T - 1, 2*genus));
      //L = CeilToZZ(sqrt(to_RR(2*L2)));
      L = CeilToZZ(sqrt(alpha*to_RR(L2)));
    }
    A = to_RR(E + L2) - power(T + 1,2*genus) ;
    if(A>0){
      L2 = FloorToZZ(power(T + 1, 2*genus)) - E;
      //L = CeilToZZ(sqrt(to_RR(2*L2)));
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

    // Finding |SUM(x_j^{lambda+1}|
    if(deg(D)%3)
      psi3 = 2.0;
    else
      psi3 = RR();

    SvSum = ZZ();
    // Finding |SUM_{v|(lambda+1), v != lambda +1}vS_v((lambda+1)/v)|
    for(j=1; j<i; j++){
      if(j%i == 0){
	if((i/j)%3)
	  SvSum += j*SvCache[2*(j-1)];
	else
	  SvSum += j*SvCache[2*j-1];
      }
    }
    psi3 += to_RR(abs(SvSum));
    psi3 /= (i*power(Q, i));

    psi3 += (2.0*genus/(i*power(T,i)) + (2*genus*T*Q*power(T, -n))/(n*(T-1)*(Q-1)) + 4.0*(Q/(Q-1))*(S/(S-1))*power(S, -n)/n);
      
    L32 = CeilToZZ(E2*expm1(psi3)+1.0/2.0);
    
    // I'm going with this guy to speed things up.
    L = CeilToZZ(sqrt(alpha*to_RR(L32)));

    // Previously, there was a correction if the new interval
    // went out of the bounds of the Hasse interval.
    // That approach changed E, so we won't do that anymore.
    
    cout<<"Hasse Interval: ["<<CeilToZZ(power(T - 1, 2*genus))<<", "<<FloorToZZ(power(T + 1, 2*genus))<<"].\n";
    //cout<<"New interval 1: ["<<E1-L12<<", "<<E1+L12<<"].\n";
    //cout<<"New interval 2: ["<<E-L2<<", "<<E+L2<<"].\n";
    cout<<"New interval 3: ["<<E-L32<<", "<<E+L32<<"].\n\n";
    //cout<<"Estimate of the class number E_1 = "<<E1<<".\n";
    cout<<"Estimate of h:                 E = "<<E<<".\n";
    //cout<<"Our bounds on the error:   L_1^2 = "<<L12<<".\n";
    //cout<<"Our bounds on the error:   L_2^2 = "<<L2<<".\n";
    cout<<"Upper bound on the error:      U = "<<L32<<".\n\n";
    //cout<<"Setting L = "<<L<<".\n\n";

    //}
  /*else{
    if(deg(D)%3 != 0) j=0;
    else j=2;
    C = power(Q,lambda);
    for(i=lambda + 1; i<= 2*lambda; i++){
      C*=Q;
      if(i%2 == 0)
	B+=(j + 2.0*genus*power(Q,i/2) + (2.0*Q/(Q-1))*(power(Q,i/2) - 1))/(i*C);
    }
  }
  */
}

// Computes the value S_v(a) = SUM_{deg(p)=v}(z1(p)^a + z2(p)^a)
// in the case that q = 1 (mod 3).
// This function can probably be improved by applying sieving
// techniques to automatically eliminate a number
// of reducible polynomials, so that expensive
// irreducibility tests are done on much fewer polynomials.

// This may also be sped up by doing random samplings of 
// irreducible polynomials and extending the results.

// Array structure: 
// [S_1(1), S_1(3), S_2(1), S_2(3), ...]

void Sv1(int v, int a){
  int loc = 2*(v-1) + (a%3==0), i, divlim;
  ZZ L, j, c;
  ZZ_pX P;

  // In this case, z1(P)^a + z2(P)^a = 2 if P!|GH and 0 if P|GH.
  if(a%3 == 0){
    divlim = v/4;
    L = ZZ();
    for(i=1; i<=divlim; i++){
      if(v%i == 0)
	L+=mu(v/i)*power(q,i);
    }
    if(v%3==0)
      L-=power(q,v/3);
    if(v%2==0)
      L-=power(q,v/2);
    L+=power(q,v);

    L/=v;

    // Subtract off the number of irreducible factors of GH of degree v.
    L-=factordegs[v-1];
    
    SvCache[loc] = 2*L;
    return;
  }

  // Otherwise,  z1(P)^a + z2(P)^a = 2 if (D|P)_3 = 1 and 
  // z1(P)^a + z2(P)^a = -1 if (D|P)_3 != 1
  P = ZZ_pX(v, 1);
  j=ZZ();

  // Loop through every monic irreducible degree v polynomial.
  // In the case that deg(P) = v = 1, every polynomial is irred.
  if(v == 1){
    splitting[0] = factordegs[0];
    // If deg(H) = 1, then H(x) = x, so x|D, and we skip
    // to P = x+1.
    if(deg(H) == 1)
      j++;
    // If there are some ramified primes of this degree, 
    // We need to check those.
    if(factordegs[0]){
      while(j<q){
	// Get the next polynomial.
	SetCoeff(P, 0, to_ZZ_p(j));
	if(IsZero(D%P)); // This test can be sped up.
	else if(chi(P) == 1){
	  // This is the complete splitting case.
	  splitting[2]++;
	  L+=2;
	}
	else{
	  // This is the inert case.
	  splitting[1]++;
	  L--;
	}
	j++;
      }
    }
    else{
      while(j<q){
	// Get the next polynomial.
	SetCoeff(P, 0, to_ZZ_p(j));
	if(chi(P) == 1){
	  // This is the complete splitting case.
	  splitting[2]++;
	  L+=2;
	}
	else{
	 // This is the inert case.
	  splitting[1]++;
	  L--;
	}
	j++;
      }
    }
  }
  // We use a nice trick to run through every irreducible 
  // degree 2 polynomial.
  else if (v == 2){
    splitting[0] = factordegs[1];
    //P = ZZ_pX(v, 1);
    // If there are some ramified primes of this degree, 
    // We need to check those.
    if(factordegs[1]){
      for(j = to_ZZ(1); j < q; j++){
	// If x^2 - j is irreducible, cycle through each (x+c)^2 - j.
	if(Jacobi(j, q) == -1){
	  SetCoeff(P, 0, to_ZZ_p(-j));
	  for(c = ZZ(); c < q; c++){
	    if(IsZero(D%P)); // This test can be sped up.
	    else if(chi(P) == 1){
	      // This is the complete splitting case.
	      splitting[2]++;
	      L+=2;
	    }
	    else{
	      // This is the inert case.
	      splitting[1]++;
	      L--;
	    }
	    // Get the next polynomial.
	    SetCoeff(P, 0, to_ZZ_p(rep(coeff(P,0))+(2*c+1)));
	    SetCoeff(P, 1, to_ZZ_p(rep(coeff(P,1))+2));
	  }
	}
      }
    }
    else{
      for(j = to_ZZ(1); j < q; j++){
	// If x^2 - j is irreducible, cycle through each (x+c)^2 - j.
	if(Jacobi(j, q) == -1){
	  SetCoeff(P, 0, to_ZZ_p(-j));
	  for(c = ZZ(); c < q; c++){
	    if(chi(P) == 1){
	      // This is the complete splitting case.
	      splitting[2]++;
	      L+=2;
	    }
	    else{
	      // This is the inert case.
	      splitting[1]++;
	      L--;
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
    // Here we have deg(P) > 2, so we test for irreducibility.
    splitting[0] = factordegs[v-1];
    if(factordegs[v-1]){
      while(IsOne(LeadCoeff(P))){
	// Formerly if(ProbIrredTest(P, iter=1)){
	if(DetIrredTest(P)){
	  if(IsZero(D%P));
	  else if(chi(P) == 1){
	    splitting[2]++;
	    L+=2;
	  }
	  else{
	    splitting[1]++;
	    L--;
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
    else{
      while(IsOne(LeadCoeff(P))){
	
	if(DetIrredTest(P)){
	  if(chi(P) == 1){
	    splitting[2]++;
	    L+=2;
	  }
	  else{
	    splitting[1]++;
	    L--;
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
  
  SvCache[loc] = L;
  return;

}

// Computes the value S_v(a) = SUM_{deg(p)=v}(z1(p)^a + z2(p)^a)
// in the case that q^v = 2 (mod 3), that is when 
// q = 2 (mod 3) and v is odd.
  

void Sv2(int v, int a){
  int loc = 2*(v-1) + (a%2==1), i, divlim;
  ZZ L;

  //  cout<<"v = "<<v<<", n/v = "<<a;

  // If v is odd and a is odd, then z1(p)^a + z2(p)^a = 0.

  if((v%2) && (a%2)){
    SvCache[loc] = ZZ();
    // cout<<", and S = "<<ZZ()<<endl;
    return;
  }
 
  // If v is even, or v is odd and a is even, then z1(p)^a + z2(p)^a = 2.
  // So we just need to calculate the number Lq(n) of monic irreducible
  // polynomials of degree n over Fq, subtract the number of irreducible
  // factors of G and H of degree v, multiply by 2, and return.
  divlim = v/4;
  L = ZZ();
  for(i=1; i<=divlim; i++){
    if(v%i == 0)
      L+=mu(v/i)*power(q,i);
  }
  if(v%3==0)
    L-=power(q,v/3);
  if(v%2==0)
    L-=power(q,v/2);
  L+=power(q,v);

  L/=v;

  // cout<<"Number of irreducible polys of degree "<<v<<" = "<<L<<endl;
  
  // Subtract off the number of irreducible factors of GH of degree v.
  L-=factordegs[v-1];

  L*=2;
  SvCache[loc] = L;
  //  cout<<", and S = "<<L<<endl;
  return;
}

// The Moebius mu function.

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

// Computes the character [D|P]_3.
// Follows Algorithm 6.2 in Class Number Algorithms in Cubic Function Fields
// by Scheidler and Stein.
ZZ_p chi(ZZ_pX P){
  int f, p;
  ZZ_p e = to_ZZ_p(1);
  ZZ_pX F = D, T;

  //GCD(T, F, P);
  
  //cout<<P<<" "<<F<<" "<<T<<endl;

  //if(!IsOne(T))
  //  return ZZ_p();

  p = deg(P)%3;
  while(deg(F)>0){
    F%=P;
    f = deg(F)%3;
    e*=(power(rec3(LeadCoeff(F)), p)*power(inv(rec3(LeadCoeff(P))),f));
    T = F;
    F = P;
    P = T;
    p = f;
  }  
  return(e*power(rec3(LeadCoeff(F)), deg(P)%3));

}

// Computes the cubic reciprocity of a mod q.
inline ZZ_p rec3(ZZ_p a){
  return(power(a, (q-1)/3));
}

RR min(RR a, RR b){
  return (a < b ? a : b);
}

// Returns a factor of the order h.
// Works only for rank 0 curves.

ZZ smallOrder0(){
  ZZ k=to_ZZ(3), factor = to_ZZ(1); //, ordg;
  //cubic_ideal g, gs, I, z, temp;
  //int i, tests = 4;
  //ZZ numSteps = CeilToZZ(sqrt(to_RR(L)/to_RR(factor)));
  //long j, place=-1, hashsize = 3*to_long(numSteps);
  vec_pair_ZZ_pX_long factors;
  //cubic_ideal baby[hashsize];
  //ZZ baby_loc[hashsize];

  //  g.G = G;       g.H = H;
  //I.G = G;       I.H = H;
  //gs.G = G;     gs.H = H;
  //z.G = G;       z.H = H;
  //temp.G = G; temp.H = H;

  // Ramified primes have order 3. 
  // So we know immediately that 3|h.
  /*
  if(deg(H) == 1)
    factor = to_ZZ(3);
  if(deg(H) > 1){
    factors = berlekamp(inv(LeadCoeff(H))*H);
    for(i=0; i<factors.length(); i++){
      if(deg(factors[i].a) <= genus)
	factor = to_ZZ(3);
      
    } 
  }
  */
  //if(IsOne(factor)){
  factors = berlekamp(inv(LeadCoeff(G))*G);
  
  //factors = berlekamp(inv(LeadCoeff(D))*D);
  if(IsOne(H)){
    power(factor, k, factors.length()-1);
  }
  else{
    power(factor, k, factors.length());
  }

  //for(i=0; i<factors.length(); i++){
    // Ramified primes have order 3. 
    // So we know immediately that 3|h.
    //if(deg(factors[i].a) <= genus)
    //  factor = to_ZZ(3);   
  /*  
      if(deg(factors[i].a) <= genus){
      z.s = factors[i].a;
      z.s2 = z.s1 = ZZ_pX(0,1);
      z.u = z.v = z.w = ZZ_pX();
      z = z*g;
      g=z;
      temp = z^factor;
      if(temp == I){
	temp = temp^k;
	if(temp == I)
	  factor*=k;
      }
    }
  } 
  */
  //cout<<"Number of factors = "<<factors.length()<<endl;

    /*
      if(factors.length() == 1){
      z.s = G;
      z.s2 = z.s1 = ZZ_pX(0,1);
      z.u = z.v = z.w = ZZ_pX();
      reduce(g, z, D);
      ordg = makeBabySteps(g, gs, baby, baby_loc, hashsize, numSteps, factor,0);
      if(!IsZero(ordg))
      factor = ordg;
      }
    
    */

  // We'll pick a few random ideals and see if they have small order.
  /*
    for(i=0; i<tests; i++){
    random(g);
    ordg = makeBabySteps(g, gs, baby, baby_loc, hashsize, numSteps, factor, 0);
    
    //cout<<"g = \n";
    //g.print();

    //temp = g^factor;
    //cout<<"g^factor = \n";
    //temp.print();

    //cout<<"factor = "<<factor<<endl;
    //cout<<"ordg = "<<ordg<<endl;
 
    // Do Giant Steps
    if(IsZero(ordg)){
    z = gs;
    k = ZZ();
    place = -1;
    while((k<numSteps)&&(place<0)){
    place = search(z, baby, hashsize);      
    temp = cubic_ideal(z);
    z = temp*gs;
    k++;
    }
    // If we've found a match.
    if(place >= 0){
    //cout<<factor<<" "<<k<<" "<<place<<endl;
    factor*=(numSteps*k+ baby_loc[place]);
    //temp = g^factor;
    //cout<<"g^"<<factor<<" = \n";
    //temp.print();
    }
    }
    else{
    //cout<<factor<<" "<<ordg<<endl;
    factor = ordg;
    }
    // Reset the baby step hash table.
    for(j=0; j<hashsize; j++)
    baby[j] = I;
    }
  */
  return factor;

}

// Checks to see if ord = 1 mod 3 is a possible order
// returns 1 if it is a possible order
// returns 0 if it is not.
int checkPrimes(ZZ &ord){
  ZZ test = ord;
  //ZZ extra=to_ZZ(1);
  int k=0;
  int i;
 

  while(test%3 == 0){
    test/=3;
  }

  if(test%3 == 2)
    return 0;

  // Check all the primes up to 200.
  for(i=0; (i<45) && !IsOne(test); i++){
    k=0;
    while(test%primes[i] == 0){
      test/=primes[i];
      k++;
    }
    if(primes[i]%3 == 2){
      if(k%2)
	//extra*=primes[i];
	return 0;
    }

    //if(test%3 == 1){
    //  factor*=extra;
    //  return;
    //}
  }

  if(test%3 == 1){
    //factor*=extra;
    return 1;
  }

  return 0;

  //if(ProbPrime(test)){
  //  extra*=test;
  //  factor*=extra;
  //  return;
  //}

  // If it gets here, then factor has a couple factors larger than 200,
  // an odd number of which are = 2 mod 3. We should improve this to pull
  // out more factors
    
}

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


// Baby-step, Giant-step algorithm for unit rank 0.
// Input: factor - an integer that we know is a factor of the order h.
// 
// Knowing so speeds up the algorithm by a factor of sqrt(so).

void bsgs0(ZZ& factor){
  int found = 0, genTrials=0, test1=200, test2=10, i;
  long hashsize, place1, place2;
  cubic_ideal g; // A subgroup generator.
  cubic_ideal z, temp, last, zinv;
  cubic_ideal gs; // The giant step.
  cubic_ideal I;  // An identity element.
  ZZ start, end = E+L2; // = E-L2, end=E+L2;
  ZZ order, order2, giantstep;
  ZZ ordg;
  ZZ numSteps;// = CeilToZZ(to_RR(L)/sqrt(to_RR(factor)));
  ZZ stepLen = factor;
  ZZ k = ZZ(), percent, r;
  //RR A, T = sqrt(to_RR(q));
  time_t set, bstime, gstime;
  int pComplete = 0;
  double n;   // n = #time for one giant step/#time for one baby step.

  if(genus == 3){
    n = 1.38;
  }
  else if(genus == 4){
    n = 1.51;
  }
  else if( genus == 5){
    n = 1.393;
  }
  else if( genus == 6 ){
    n = 1.615;
  }
  else if( genus == 7 ){
    n = 1.54;
  }
  else{
    n = 1.615; 
  }

  //cubic_ideal *baby;
  //ZZ *baby_loc;
  //ZZVec baby_loc;

  // Testing the whole Hasse Interval if needed.
  //start = to_ZZ(430443912);
  //end = to_ZZ(656135513);
  //numSteps = CeilToZZ(sqrt(to_RR(end-start)/to_RR(factor)));

  order0 = ZZ();
  order1 = ZZ();

  z.G = G;
  z.H = H;
  temp.G = G;
  temp.H = H;
  gs.G = G;
  gs.H = H;

  // The number of Baby Steps. It will be much nicer if they are even.
  numSteps = CeilToZZ(to_RR(L)*sqrt(n)/sqrt(to_RR(stepLen)));
  numSteps += (numSteps%2);
  //cout<<"numSteps = "<<numSteps<<endl;

  hashsize = NextPrime(4*to_long(numSteps));
  hashbits = NumBits(hashsize);

  //cubic_ideal baby[hashsize];
  
  // We're only going to store a couple coefficients of s 
  // from the ideal to save time and space.

  //cout<<"hashsize = "<<hashsize<<" "<<q.size()<<" "<<L.size()<<endl;
  ZZVec baby1 = ZZVec(hashsize, q.size());  // store s(0)
  ZZVec baby2 = ZZVec(hashsize, q.size());  // store s(1)
  ZZVec baby_loc = ZZVec(hashsize, L.size());

  if((&baby1 == NULL) || (&baby2 == NULL) || (&baby_loc == NULL) ){
    cout<<"Not enough memory available to create baby steps.\n";
    cout<<"Hash table size = "<<hashsize<<".\n"<<"Total memory required = "<<3*hashsize*sizeof(ZZ)<<" bytes.\n\n";
    return;
  }

  g.G = G;
  g.H = H;
  
  //D = G*sqr(H);
 
  /**************/
  /* BABY STEPS */
  /**************/

  cout<<"Phase 3: Making Baby Steps ... \n";
  set = time(NULL);
  do {
    // Initialize the baby step hash table.
    //baby = (cubic_ideal *)malloc(hashsize*sizeof(cubic_ideal));
    //baby_loc = ZZVec(hashsize, NumBits(numSteps));
    //baby_loc = (ZZ *)malloc(hashsize*sizeof(ZZ));

    //cout<<baby<<endl;
    //cout<<baby_loc<<endl;

    for(long i=0; i<hashsize; i++){
      baby_loc[i] = N1;
      //baby[i] = I;
      clear(baby1[i]);
      clear(baby2[i]);
    }

    // Get a random (sub)group generator.

    random(g);

    //cout<<"numSteps = "<<numSteps<<endl;

    //g.s = g.s1 = ZZ_pX(3,0)+ZZ_pX(2,1)+ZZ_pX(1,15)+ZZ_pX(0,16);
    //g.s2 = ZZ_pX(1,0)+ZZ_pX(0,1);
    //g.u = ZZ_pX();//ZZ_pX(1,13)+ZZ_pX(0,3);
    //g.v = ZZ_pX(2,0)+ZZ_pX(1,6)+ZZ_pX(0,12);
    //g.w = ZZ_pX(1,8)+ZZ_pX(0,8);
    
    g.print();

    ordg = makeBabySteps(g, baby1, baby2, baby_loc, hashsize, numSteps, stepLen, 1);

    //cout<<endl;
    //for(long i=0; i<hashsize; i++)
    //  cout<<baby_loc[i]<<" ";
    //cout<<endl;
    //cout<<"max Collision = "<<maxCollision<<endl;

    // Check other generators to see if this is the order.
    // Or at least try to find factors of the order.
    if(!IsZero(ordg)){
      genTrials = 0;
      while(genTrials < test2){
	random(z);
	temp = z^order;
	if(temp != I)
	  break;
	else
	  genTrials++;
      }
    
      // We probably got the jackpot here.
      if(genTrials == test2){
	order = ordg;
	goto end;
      }
      else{
	// This could give us a factor of g.
	if(!IsZero(order0) && (!IsZero(order1))){
	  
	  ordg = order0-order1;
	  temp = g^ordg;
	  // We found a factor.
	  if(IsOne(temp.s)){
	    factor = ordg;
	
	    stepLen = factor;
	    numSteps = CeilToZZ(to_RR(L)/sqrt(to_RR(stepLen)));
	    numSteps += (numSteps%2);
	    ordg = ZZ();
	  }
	}
      }
    }
    else{
      // We could have a factor of g.
      if(!IsZero(order0) && (!IsZero(order1))){
	
	ordg = order0-order1;
	temp = g^ordg;
	// We found a factor.
	if(IsOne(temp.s)){
	  factor = ordg;
	  
	  stepLen = factor;
	  numSteps = CeilToZZ(to_RR(L)/sqrt(to_RR(stepLen)));
	  numSteps += (numSteps%2);
	ordg = ZZ();
	}
      }
    }
  } while(!IsZero(ordg));
  
  genTrials = 0;
  bstime = time(NULL);
  cout<<"\n"<<"Baby Steps completed in "<<bstime-set<<" seconds. \nMaking Giant Steps ... \n";

  /***************/
  /* GIANT STEPS */
  /***************/

  // Set the giant step distance and the first giant step.
  cout<<stepLen<<" "<<numSteps<<endl;
  if(irred && (deg(H)==0)){
    giantstep = stepLen*numSteps;
    gs = g^giantstep;
  }
  else{
    giantstep = stepLen*numSteps/2;
  }
  gs = g^giantstep;
  z = cubic_ideal(gs);

  //cout<<"giantstep = "<<giantstep<<endl;
  
  percent = L2/(100*giantstep);
  if(percent==0) percent=1;
  pComplete = 0;

  do{
    do{
      place1 = search(z, baby1, baby2, baby_loc, hashsize);

      if(place1 >= 0)
	found = 1;
      
      // Search the hash table for the inverse of z.
      //cout<<"z = "<<z.is_valid()<<endl;
      //z.print();
      
      inverse(temp, z);
      //cout<<"temp = "<<temp.is_valid()<<endl;
      //temp.print();
      reduce(zinv, temp, D);
      
      place2 = search(zinv, baby1, baby2, baby_loc, hashsize);
      
      if(place2 >= 0)
	found = 2;
      
      if(IsZero(k%percent)){
	printf("\r");
	printf("Giant steps are %d %% complete. ", pComplete);

	cout<<flush;
	pComplete++;
      }

      temp = cubic_ideal(z);
      //cout<<"temp = "<<temp.is_valid()<<endl;
      //temp.print();
      //cout<<"gs = "<<gs.is_valid()<<endl;
      //gs.print();

      z = temp*gs;
      k++;

      //cout<<"I am here! "<<k<<endl;
    }while((place1 < 0) && (place2 < 0));    

    //cout<<"I am here!"<<endl;
    //if(k==numSteps){
    //cout<<"Went overboard!\n\n";
    //cout<<"Up to: "<<start + factor*numSteps*(k-1)<<endl;
    //gstime = time(NULL);
    //cout<<"Giant steps \"completed\" in "<<gstime-bstime<<" seconds.\n\n";
    //return;
    // }
    //cout<<"k = "<<k<<" places = "<<place1<<" "<<place2<<endl;
    // At this point we've found a match.
    //start = E - (E%factor);
    start = E - (E%stepLen);

    //if(irred && (deg(H) == 0)){
    //  start++;
    //}

    //cout<<"start = "<<start<<" giantstep = "<<giantstep<<" k = "<<k<<" factor = "<<factor<<" location = "<<baby_loc[place]<<endl;
    if(place1 >= 0){
      //order2 = start + giantstep*(k-1);
      order = start - giantstep*k + baby_loc[place1];
    }
    else{
      //order2 = start + giantstep*(k-1);
      order = start + giantstep*k + baby_loc[place2];
    }
    //cout<<"k = "<<k<<endl;
    cout<<"Match found.\n"<<"Order = "<<order<<endl;

    //cout<<place1<<" "<<place2<<" z.s = "<<temp.s<<endl;

    // Check to see if it is valid.
    if(checkPrimes(order)){
    
      // Test several elements to see if this is the order.
      genTrials = 0;
      while(genTrials < test2){
	//cout<<"Getting random ideal: "<<endl;
	random(g);
	//g.s = g.s1 = ZZ_pX(3,0)+ZZ_pX(2,1)+ZZ_pX(1,15)+ZZ_pX(0,16);
	//g.s2 = ZZ_pX(1,0)+ZZ_pX(0,1);
	//g.u = ZZ_pX();//ZZ_pX(1,13)+ZZ_pX(0,3);
	//g.v = ZZ_pX(2,0)+ZZ_pX(1,4)+ZZ_pX(0,8);
	//g.w = ZZ_pX(1,12)+ZZ_pX(0,12);

	//g.s = ZZ_pX(0, 10)+ZZ_pX(1, 6)+ZZ_pX(2, 17)+ZZ_pX(3, 7)+ZZ_pX(4, 1);
	//g.s1 = g.s2 = ZZ_pX(0,1);
	//g.u = ZZ_pX(0,3)+ZZ_pX(1, 16);
	//g.v = ZZ_pX(0,12)+ZZ_pX(1, 12)+ZZ_pX(2, 12)+ZZ_pX(3, 18);
	//g.w = ZZ_pX();
	//cout<<"Valid ideal? "<<g.is_valid()<<endl;
	//g.print();
	temp = g^order;
	//cout<<"Is g^h OK? "<<temp.is_valid()<<endl;
	//temp.print();
	//ZZ_pX D = G*sqr(H);
	//reduce(g, temp, D);
	//cout<<"reduced:"<<endl;
	//g.print();
	// This isn't the order. Keep searching.
	if(temp != I){
	  cout<<"Found the wrong order.\n\n";
	  found = 0;
	  break;
	}
	else
	  genTrials++;
      }
    }
    else{
      cout<<"Found the wrong order.\n\n";
      found = 0;
    }

    // If two matches were found, and the first was bad 
    // then check the other match.
    if((place1 >= 0) && (place2 >= 0) && (!found)){
      found = 2;
      order = start + giantstep*k + baby_loc[place2];
      if(checkPrimes(order)){
	genTrials = 0;
	while(genTrials < test2){
	  random(g);
	  temp = g^order;
	  // This isn't the order. Keep searching.
	  if(temp != I){
	    cout<<"Found the wrong order.\n\n";
	    found = 0;
	    break;
	  }
	  else
	    genTrials++;
	}
      }
      else{
	cout<<"Found the wrong order.\n\n";
	found = 0;
      }
    }
  } while(!found);

  gstime = time(NULL);
  cout<<endl;
  cout<<"Giant steps completed in "<<gstime-bstime<<" seconds.\n\n";

 end:

  cout<<"The order is h = "<<order<<".\n\n";
  cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<"    |h-E|   = "<<abs(order-E)<<endl; 
  cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  cout<<"|h-E2|/L3^2 = "<<to_RR(abs(order-E))/to_RR(L32)<<endl;
  
  // Checking.
  //temp = cubic_ideal();
  //z = g^order;
  //if(z != temp){
  //  cout<<"Error: The order h is incorrect.\n\n";
  //}

  return;
}

// Input: g: A (sub)group generator.
//        baby1,2: A hash table to store baby steps.
//        loc: An array that gives the power of g stored in the array.
//        hashsize: the size of the hashtable.
//        steps: number of baby steps to compute.
//        stepLen: The length of a baby step.
// Return: either the order of g, or 0 if the identity was not encountered. 

ZZ makeBabySteps(cubic_ideal &g, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize, ZZ& steps, ZZ& stepLen, int verbose){
  // The identity and a random (sub)group generator.
  cubic_ideal h, z, step;
  long j, index;
  ZZ i = to_ZZ(1), start; 
  ZZ steps0;
  ZZ percent = steps/100;
  if(percent==0)
    percent++;
  int pComplete = 0;

  h.set_G(G);
  h.set_H(H);
  z.set_G(G);
  z.set_H(H);
  step.G = G;
  step.H = H;

  // We'll be working in equivalence classes mod 3.
  stepLen*=3;

  start = E - (E%stepLen);

  // Make the baby steps.
  // Begin at g^E and go to g^{E+steps-1}
  //cout<<"Finding g^start.\n";
  //cout<<"Is this a good ideal? "<<g.is_valid()<<endl;
  //g.print();
  h = g^start;

  //cout<<"Is this a good ideal? "<<h.is_valid()<<endl;
  //h.print();

  step = g^stepLen;

  //cout<<"Is this a good ideal? "<<step.is_valid()<<endl;
  //step.print();

  // If D(x) is irreducible, h = 1 (mod 3). Work in this equivalence class.
  if(irred && (deg(H) == 0)){
    start++;

    z = h;
    h = z*g;   

    // Work in the equivalence classes 1 mod stepLen.
    for(i=0; i<steps; i++){
      
      if(verbose){
	if(IsZero(i%percent)){
	  printf("\r");
	  printf("Baby steps are %d %% complete.", pComplete);
	  cout<<flush;
	  pComplete++;
	}
      }

      // If the baby step is the identity, then start + i*stepLen + 1
      // is a multiple of the order of g, but not necessarily the
      // group order.
      // Check small primes, especially those = 2 mod 3.
      if(IsOne(h.s)){
	if(!IsZero(order0))
	  order1 = order0;
	order0 = start+i*stepLen + 1;
	
	cout<<"order = "<<order0<<endl;
	if(checkPrimes(order0)){
	  cout<<"Passed checkPrimes\n";
	  return (order0);
	}
      }

      //cout<<h.s<<endl;
      
      // Otherwise, we insert the value into the hash table.
      insert(h, baby1, baby2, loc, stepLen*i+1, hashsize); 

      // Make the next baby step.
      //h.print();

      z = cubic_ideal(h);

      h = z*step;

      //if(!h.is_valid()){
      //z.print();
      //step.print();
      //h.print();
      //}
    }
  }
  else{
    ZZ stepLenDiv3 = stepLen/3;
    
    // h covers 0 mod stepLen
    // h1 covers 1 mod stepLen
    cubic_ideal h1;
    h1.set_G(G);
    h1.set_H(H);
    
    z = g^stepLenDiv3;
 
    h1 = h*z;    

    // steps is always even. We begin on a = 0 (mod stepLen)
    // location in the nonirreducible case.
    steps0 = steps/2;

    percent/=2;
    if(percent==0) percent++;

    //cout<<"start = "<<start<<" "<<start%3<<" "<<stepLen<<endl;
    //cout<<steps0<<" "<<steps1<<" "<<steps<<endl;
  
    for(i=0; i<steps0; i++){
    
      if(verbose){
	if(IsZero(i%percent)){
	  printf("\r");
	  printf("Baby steps are %d %% complete.", pComplete);
	  cout<<flush;
	  pComplete++;
	}
      }
      
      // If we've nailed it, return the order.
      if(IsOne(h.s)){
	if(!IsZero(order0)){
	  order1 = order0;
	  order0 = start+i*stepLen;
	  if(checkPrimes(order0)){
	    return (order0);
	  }
	}
      	else{
	  order0 = start+i*stepLen;
	  if(checkPrimes(order0)){
	    return (order0);
	  }
	}
      }
      
      //cout<<h.s<<endl;

      //cout<<"Is h a good ideal? "<<h.is_valid()<<endl;
      //h.print();
      //cout<<"Is h1 a good ideal? "<<h1.is_valid()<<endl;
      //h1.print();

      // Otherwise, we insert the value into the hash table.
      insert(h, baby1, baby2, loc, stepLen*i, hashsize); 
      
      // Make the next baby step.
      z = cubic_ideal(h);

      //cout<<"h = "<<z.is_valid()<<endl;
      //z.print();
      //cout<<"step = "<<step.is_valid()<<endl;
      //step.print();

      h = z*step;

      // Repeat for the = 1 mod stepLen case.

      // If we've nailed it, return the order.
      if(IsOne(h1.s)){
	if(!IsZero(order0)){
	  order1 = order0;
	  order0 = start+i*stepLen+stepLenDiv3;
	  if(checkPrimes(order0)){
	    return (order0);
	  }
	}
      	else{
	  order0 = start+i*stepLen+stepLenDiv3;
	  if(checkPrimes(order0)){
	    return (order0);
	  }
	}
      }

      // Otherwise, we insert the value into the hash table.
      insert(h1, baby1, baby2, loc, stepLen*i+ stepLenDiv3, hashsize); 
      
      // Make the next baby step.
      z = cubic_ideal(h1);
      
      //cout<<"h1 = "<<z.is_valid()<<endl;
      //z.print();
      //cout<<"step = "<<step.is_valid()<<endl;
      //step.print();
            
      h1 = z*step;

    }
  }
  //cout<<"Max Collision = "<<maxCollision<<endl;
  return ZZ();
}

// The Kangaroo Method

void kangaroo0(ZZ& factor){
  int i, found = 0;
  int distbits; // The number of zero bits to determine 
                // if an ideal is distinguished.
  int distbits2;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  ZZ avgjump;  // The average of the jumps.
  ZZ randjump; // The upper bound on the random jump distance.
  ZZ start;
  cubic_ideal g, tame, wild;
  cubic_ideal temp;
  cubic_ideal jumps[50];
  ZZ jumpdistance[50];
  ZZ jumpsum = ZZ();
  
  ZZ tamedist, wilddist = ZZ(), matchdist = ZZ();

  ZZ jumptotal = ZZ(); // The number of jumps for each kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  long hashsize;
  ZZ order; // The order of the group.

  // Initialize the cubic_ideals.
  g.G = G; g.H = H;
  tame.G = G; tame.H = H;
  wild.G = G; wild.H = H;
  temp.G = G; temp.H = H;
  for(i=0; i<50; i++){
    jumps[i].G = G; jumps[i].H = H;
  }
  
  // Set the average jump distance, the jump distances, 
  // and the starting position for tame.
  if(IsOne(factor)){
    // In this case h = 1 mod (3), so the jump distances must be
    // multiples of 3. 
    avgjump = 3*L;
    
    //avgjump-=(avgjump%3);
    randjump = 2*avgjump/3;
    jumpdistance[49]=0;
    for(i=0; i<49; i++){
      jumpdistance[i] = 3*(RandomBnd(randjump)+1);
      
      jumpsum += jumpdistance[i];
      while( (jumpsum < (i-1)*avgjump) || (jumpsum > (i+1)*avgjump) ){
	jumpsum -= jumpdistance[i];
	jumpdistance[i] = 3*(RandomBnd(randjump)+1);
	jumpsum += jumpdistance[i];
      }
      
      jumpdistance[49] -= jumpdistance[i];
      //cout<<jumpdistance[i]<<" ";
    }
    //cout<<jumpdistance[49]<<endl;
    jumpdistance[49] += 50*avgjump;
    //cout<<jumpdistance[49]<<endl;

    start = E - E%3 + 1;
  }
  else{
    RoundToZZ(avgjump, to_RR(L)*sqrt(1.5)*sqrt(to_RR(factor)));
    avgjump-=(avgjump%factor);
    randjump = 2*avgjump/factor;
    jumpdistance[49]=0;
    for(i=0; i<49; i++){
      jumpdistance[i] = factor*(RandomBnd(randjump)+1);
      jumpsum += jumpdistance[i];
      while( (jumpsum < (i-1)*avgjump) || (jumpsum > (i+1)*avgjump) ){
	jumpsum -= jumpdistance[i];
	jumpdistance[i] = factor*(RandomBnd(randjump)+1);
	jumpsum += jumpdistance[i];
      }
      jumpdistance[49] = -jumpdistance[i];
    }
    
    jumpdistance[49] += 50*avgjump;

    start = E - E%factor;
  }

  // Set the number of bits for a distinguished ideal.
  distbits = (NumBits(avgjump)-1)/2;

  distbits += (distbits%2);
  distbits2 = distbits/2;

  // Set the hashtables and their size.
  hashsize = NextPrime(to_long(10*RightShift(avgjump, distbits)));
  hashbits = NumBits(hashsize);

  //cout<<avgjump<<" "<<distbits<<" "<<hashsize<<endl;

  // The hash tables for the tame kangaroo.

  ZZVec dpt1 = ZZVec(hashsize, q.size());  // store s(0)
  ZZVec dpt2 = ZZVec(hashsize, q.size());  // store s(1)
  ZZVec dptd = ZZVec(hashsize, E.size());  // store the distance

  // The hash tables for the wild kangaroo.

  ZZVec dpw1 = ZZVec(hashsize, q.size());  // store s(0)
  ZZVec dpw2 = ZZVec(hashsize, q.size());  // store s(1)
  ZZVec dpwd = ZZVec(hashsize, E.size());  // store the distance

  // Initialize the hashtables.
  for(long j=0; j<hashsize; j++){
    dptd[j] = dpwd[j] = N1;
  }

  // Set the jumps. First generate a random generator.
  random(g);

  for(i=0; i<50; i++){
    jumps[i] = g^jumpdistance[i];
  }

  // Set the tame kangaroo. (The wild kangaroo is initialized at I.)
  tame = g^start;
  tamedist = start;
  cout<<"BOING! BOING! BOING! BOING! BOING! BOING! BOING! BOING! BOING! \n\n";

  while(!found){

    // Make the jumps.
    i = vmap(tame);
    temp = cubic_ideal(tame);
    tame = temp*jumps[i];
    tamedist += jumpdistance[i];

    i = vmap(wild);
    temp = cubic_ideal(wild);
    wild = wild*jumps[i];
    wilddist += jumpdistance[i];

    /***********************************/
    /* Check for distinguished points. */
    /***********************************/

    /***************************/
    /* Check the tame kangaroo */
    /***************************/

    roocoeff = rep(coeff(tame.s,0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(tame.s,1));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(tame, 0, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, tamedist, matchdist, hashsize); 
	distpoints++;
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
	
	// WOO HOO! There's a match!
	if(found){
	  order = abs(tamedist-matchdist);
	  // Check to make sure the order is correct.
	  // If not, insert the ideal.
	  if(!check(order)){
	    found = 0;
	    rooinsert(tame, dpt1, dpt2, dptd, tamedist, hashsize);
	  }
	}
      }
    }

    /***************************/
    /* Check the wild kangaroo */
    /***************************/

    roocoeff = rep(coeff(wild.s,0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
       roocoeff = rep(coeff(wild.s,1));
       trunc(dptest, roocoeff, distbits2);
       if(IsZero(dptest)){
	 // Check for a match, or insert otherwise.
	 found = distpoint(wild, 1, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, wilddist, matchdist, hashsize);
	 distpoints++;
	 printf("\r");
	 printf("Distinguished points: %ld", distpoints);
	 cout<<flush;
	 
	 // WOO HOO! There's a match!
	 if(found){
	   order = abs(matchdist-wilddist);
	   // Check to make sure the order is correct.
	   // If not, insert the ideal.
	   if(!check(order)){
	     found = 0;
	     rooinsert(wild, dpw1, dpw2, dpwd, wilddist, hashsize);
	   }
	 }
       }
    }
    jumptotal+=2;
  }
  cout<<"\n\n";
  cout<<"The order is h = "<<order<<".\n\n";
  cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<"    |h-E|   = "<<abs(order-E)<<endl; 
  cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  cout<<"|h-E2|/L3^2 = "<<to_RR(abs(order-E))/to_RR(L32)<<endl;
  if(IsOne(factor)){
    cout<<"Expected kangaroo jumps:       "<<4*L<<endl;
    cout<<"Expected distinguished points: "<<RightShift(L, distbits-2)<<endl;
  }
  else{
    cout<<"Expected kangaroo jumps:       "<<4*avgjump/factor<<endl;
    cout<<"Expected distinguished points: "<<RightShift(avgjump/factor, distbits-2)<<endl;
  }
  cout<<"Total kangaroo jumps:          "<<jumptotal<<endl;
  cout<<"Total distinguished points:    "<<distpoints<<endl;
}

// A function mapping a kangaroo to {0, 1, ..., 49}

inline int vmap(cubic_ideal &A){
  return(rep(coeff(A.u,0))%50);
}

// Searches the distinguished point arrays for a match, and inserts otherwise.
// torw = 0 if roo is tame
// torw = 1 if roo is wild
// Return 0 if there is no match, and 1 if there is.
int distpoint(cubic_ideal &roo, int torw, ZZVec &dpt1, ZZVec &dpt2, ZZVec &dptd, ZZVec &dpw1, ZZVec &dpw2, ZZVec &dpwd, ZZ distance, ZZ &match, long size){
  long index;

  // If the roo is tame, search for a match in the wild dp tables.
  if(torw == 0){
    index = roosearch(roo, dpw1, dpw2, dpwd, size);
    // If we didn't find a match, insert the tame kangaroo into the hash table.
    if(index == -1){
      rooinsert(roo, dpt1, dpt2, dptd, distance, size);
      return 0;
    }
    // Sweet! We found a match!
    else{
      match = dpwd[index];
      return 1;
    }
  }

  // If the roo is wild, search for a match in the tame dp tables.
  if(torw == 1){
    index = roosearch(roo, dpt1, dpt2, dptd, size);
    // If we didn't find a match, insert the wild kangaroo into the hash table.
    if(index == -1){
      rooinsert(roo, dpw1, dpw2, dpwd, distance, size);
      return 0;
    }
    // Sweet! We found a match!
    else{
      match = dptd[index];
      return 1;
    }
  }

}

// The hash function for kangaroos.

inline long roohash(cubic_ideal &A, long size){
  return( ((rep(coeff(A.v,0))<<hashbits) + rep(coeff(A.s,1)))%size );
}

// Insert a kangaroo into its hash table.

inline void rooinsert(cubic_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, ZZ distance, long size){
  long index, j;

  index = roohash(A, size);
      
  // Collision resolution - quadratic probing.
  //while(baby[index]!=I){
  j = 1;
  while(loc[index] != N1){
    index+=(j*j);
    index%=size;
    j++;
  }
  
  if(j>maxCollision)
    maxCollision = j;
  
  // Insert the kangaroo.
  roo1[index] = rep(coeff(A.s,1));
  roo2[index] = rep(coeff(A.s,2));
  loc[index] = distance;
}

// Returns the position of A in the hash table, or -1 if it is not found.

inline long roosearch(cubic_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, long hashsize){
  long j, index;
  int i;

  index = roohash(A, hashsize);

  if(loc[index] == N1)  
    return(-1);
  
  if(roo1[index] == rep(coeff(A.s,1))){
    if(roo2[index] == rep(coeff(A.s,2)))
      return index;
  }
  
  // Collision resolution - quadratic probing.
  j = 1;
  for(i=0; i<maxCollision; i++){
    index+=(j*j);
    index%=hashsize;
    j++;
    if(loc[index] == N1)
      return(-1);
    if(roo1[index] == rep(coeff(A.s,1))){
      if(roo2[index] == rep(coeff(A.s,2)))
	return index;
    }
  }
  return(-1);
}


// The hash function.

inline long hash(cubic_ideal &A, long size){
  //long val = 1;
  long value = 0;
  int i, d = deg(A.u);

  if(IsZero(A.u)){
    if(IsZero(A.v)){
      if(IsOne(A.s))
	return 1;
      //return val;
      else{
	return( ((rep(coeff(A.s,0))<<hashbits) + rep(coeff(A.s,1)))%size );
      }
    }
    else{
      d = deg(A.v);
      value = ((rep(coeff(A.v, 0)) << 2) + rep(coeff(A.s,0))) % size;
      for(i=1; i<=d; i++){
	value = ((value << (hashbits/2)) + rep(coeff(A.v, i))) % size; 
      }
      return(value);
      //return( ((rep(coeff(A.v,0))<<hashbits) + rep(coeff(A.s,0)))%size );
    }
  }
  else{ 
    value = ((rep(coeff(A.u, 0)) << 2) + rep(coeff(A.v,0))) % size;
    for(i=1; i<=d; i++){
      value = ((value << (hashbits/2)) + rep(coeff(A.u, i))) % size; 
    }
    return(value);
    //return( ((rep(coeff(A.u,0))<<hashbits) + (rep(coeff(A.v,0))<<(hashbits/2)) + rep(coeff(A.s,1)) )%size );
  }

  //return val;
}

// The hash function for infrastructure ideals.

inline long hash(infrastructure_ideal &A, long size){
  //long val = 1;
  long value=0;
  int i, d = deg(A.mu0);
  
  //cout<<hashbits<<" "<<d<<": ";
  for(i=0; i<=d; i++){
    value = ((value << (hashbits)) + rep(coeff(A.mu0, i))) % size; 
    //cout<<value<<" ";
  }
  d = deg(A.mu1);
  for(i=0; i<=d; i++){
    value = ((value << (hashbits)) + rep(coeff(A.mu1, i))) % size; 
    //cout<<value<<" ";
  }
  //cout<<value<<" "<<size<<endl;
  return(value);
}

// Insert a cubic ideal into the baby step hash table.

inline void insert(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, ZZ val, long size){
  long index, j;

  // Otherwise, we insert the value into the hash table.
  index = hash(A, size);
  
  // Collision resolution - quadratic probing.
  //while(baby[index]!=I){
  j = 1;
  while(loc[index] != N1){
    index+=(j*j);
    index%=size;
    j++;
  }
  
  if(j>maxCollision)
    maxCollision = j;
  
  // Insert the ideal.
  baby1[index] = baby1val(A);
  baby2[index] = baby2val(A);
  
  loc[index] = val;
}

// Insert an infrastructure ideal into the baby step hash table.

inline void insert(infrastructure_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long size){
  long index, j;

  // Otherwise, we insert the value into the hash table.
  index = hash(A, size);
  //cout<<index<<" ";
  // Collision resolution - quadratic probing.
  //while(baby[index]!=I){
  j = 1;
  while(loc[index] != N1){
    index+=(j*j);
    index%=size;
    j++;
  }
  //cout<<index<<" ";
  if(j>maxCollision)
    maxCollision = j;
  
  // Insert the ideal.
  baby1[index] = baby1val(A);
  baby2[index] = baby2val(A);
  //cout<<"vals = "<<baby1[index]<<" "<<baby2[index]<<endl;
  loc[index] = A.d0;
}

// Returns the position of A in the hash table, or -1 if it is not found.

inline long search(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize){
  long j, index;
  int i;
  ZZ bvalue1, bvalue2;

  index = hash(A, hashsize);
  bvalue1 = baby1val(A);
  bvalue2 = baby2val(A);
  
  if(loc[index] == N1)  
    return(-1);

  if(baby1[index] == bvalue1){ //rep(coeff(A.s,0))){
    if(IsOne(A.s)){
      if(IsZero(baby2[index]))
	return index;
    }
    else{
      if(baby2[index] == bvalue2) //rep(coeff(A.s,1)))
	return index;
    }
  }
  
  // Collision resolution - quadratic probing.
  j = 1;
  for(i=0; i<maxCollision; i++){
    //while(baby[index]!=I){
    //while(!IsZero(baby[index].s)){
    index+=(j*j);
    index%=hashsize;
    j++;
    //if(baby[index]==I)
    if(loc[index] == N1)
      return(-1);
    //if(baby[index] == A)
    if(baby1[index] == bvalue1){ //rep(coeff(A.s,0))){
      if(IsOne(A.s)){
	if(IsZero(baby2[index]))
	  return index;
      }
      else{
	if(baby2[index] == bvalue2)  //rep(coeff(A.s,1)))
	  return index;
      }
    }
  }
  
  return(-1);
}

// Returns the position of A in the hash table, or -1 if it is not found.

inline long search(infrastructure_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize){
  long j, index;
  int i;
  ZZ bvalue1, bvalue2;

  index = hash(A, hashsize);
  bvalue1 = baby1val(A);
  bvalue2 = baby2val(A);
  //cout<<"index = "<<index<<" vals = "<<bvalue1<<" "<<bvalue2<<" ";
  if(loc[index] == N1)  
    return(-1);

  if(baby1[index] == bvalue1){ //rep(coeff(A.s,0))){
    //if(IsOne(A.d)){
    // if(IsZero(baby2[index]))
    //	return index;
    //}
    //else{
    if(baby2[index] == bvalue2) //rep(coeff(A.s,1)))
      return index;
      //}
  }
  
  // Collision resolution - quadratic probing.
  j = 1;
  for(i=0; i<maxCollision; i++){
    //while(baby[index]!=I){
    //while(!IsZero(baby[index].s)){
    index+=(j*j);
    index%=hashsize;
    //cout<<index<<" ";
    j++;
    //if(baby[index]==I)
    if(loc[index] == N1)
      return(-1);
    //if(baby[index] == A)
    if(baby1[index] == bvalue1){ //rep(coeff(A.s,0))){
      //if(IsOne(A.s)){
      //if(IsZero(baby2[index]))
      //  return index;
      //}
      //else{
      if(baby2[index] == bvalue2)  //rep(coeff(A.s,1)))
	return index;
      //}
    }
  }
  //cout<<"index = "<<index<<endl;
  return(-1);
}

inline ZZ baby1val(cubic_ideal &A){

  int i, d = deg(A.s);
  ZZ value;
  

  if(IsOne(A.s))
    return(-N1); 
  else{
    value = rep(coeff(A.s,0));
    for(i=2; i<d; i+=2){
      value = (value << 4) + rep(coeff(A.s, i));
    }
    return value;
  }

}

inline ZZ baby2val(cubic_ideal &A){

  int i, d = deg(A.s);
  ZZ value;
  

  if(IsOne(A.s))
    return(ZZ()); 
  else{
    value = rep(coeff(A.s,1));
    for(i=3; i<d; i+=2){
      value = (value << 4) + rep(coeff(A.s, i));
    }
    return value;
  }

}
inline ZZ baby1val(infrastructure_ideal &A){

  int i, d = deg(A.d);
  ZZ value;
  

  if(IsOne(A.d))
    return(-N1); 
  else{
    value = rep(coeff(A.d,0));
    for(i=2; i<d; i+=2){
      value = (value << 4) + rep(coeff(A.d, i));
    }
    return value;
  }

}

inline ZZ baby2val(infrastructure_ideal &A){

  int i, d = deg(A.nu0);
  ZZ value;
  

  if(IsOne(A.nu0))
    return(ZZ()); 
  else{
    value = rep(coeff(A.nu0,0));
    for(i=1; i<d; i+=2){
      value = (value << 4) + rep(coeff(A.nu0, i));
    }
    return value;
  }

}

// Input: A multiple, h0, of the S-regulator, R^S, and a lower bound on R^S.
// Output: The S-regulator, R^S.

ZZ extract(ZZ h0, ZZ bound){
  ZZ h = h0, hx = to_ZZ(1);
  ZZ p = to_ZZ(1);
  ZZ pe = to_ZZ(1);
  int i=0;
  vec_ZZ factors;
  //factors.SetLength(1);
  infrastructure_ideal A;
  A.G = G; A.H = H;

  // In the future, I'd like to replace this with an actual integer factorization routine.

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
  //cout<<h<<" "<<bound<<endl;

  for(i=0; i<factors.length(); i++){
    if(factors[i] < (h/bound)){
      pe = p = factors[i];
      //cout<<"p^e = "<<pe<<endl;
      near(2*h/pe, A);
      while(IsOne(A.d) && IsZero(h%pe)){
	pe*=p;
	//cout<<"p^e = "<<pe<<endl;
	near(2*h/pe, A);
      }
      hx*=(pe/p);
    }
  }
  return h/hx;
}

ZZ makeBabySteps1(ZZVec& baby1, ZZVec& baby2, ZZVec& baby_deg, long hashsize, infrastructure_ideal& GS){
  ZZ i;
  int pComplete = 1;
  ZZ percent = L/100;
  infrastructure_ideal A, B, C;
  C.G = B.G = A.G = G; C.H = B.H = A.H = H; 
  reduce_basis(A);

  near(2*E, A);
  C = A;
   
  //cout<<"Computing near("<<E<<")"<<endl;
  baby_step_r1(A, B);
  //cout<<"Distance check: "<<C.d0<<" <= "<<2*E<<" <= "<<B.d0<<endl;

  insert(A, baby1, baby2, baby_deg, hashsize);
  //cout<<"i = 0"<<endl;
  //A.print();

  if(IsZero(percent))
    percent++;

  for(i=1; i<L; i++){
    if(IsZero(i%percent)){
      printf("\r");
      printf("Baby steps are %d %% complete.", pComplete);
      cout<<flush;
      pComplete++;
    }
    baby_step_r1(A, B);
    A = B;
    //if(i==284){
    //  cout<<"i = "<<i<<" of "<<L<<"; A.d = "<<A.d<<" "<<A.d0<<endl;
    //  A.print();
    //}
    insert(A, baby1, baby2, baby_deg, hashsize);
    // In the off-chance that we actually make it around the cycle...
    // ... we've found the regulator!
    if(IsOne(A.d)){
      GS = infrastructure_ideal(A);
      return(A.d0/2);
    }
    if(C == A){
      GS = A;
      return((A.d0 - C.d0)/2);
    }
  }
  //cout<<"\nCompleted baby steps. A.d0 = "<<A.d0<<endl;
  //A.print();
  //GS = infrastructure_ideal(A);
  near(A.d0 - C.d0 - 2, GS);
  //cout<<"GS.d0 = "<<GS.d0<<" GS = "<<endl;
  //GS.print(); 
  return ZZ();
}


// Baby-step, Giant-step algorithm for unit rank 1.

void bsgs1(){
  ZZ R, order;
  infrastructure_ideal A, B, I, GS;
  time_t set, bstime, gstime;
  long place1, place2, k;
  int found = 0;
  ZZ_pX D = G*H*H;
  ZZ_pX s, s1, s2;
  
  double n;
  cubic_ideal T;
  
  if (genus == 3){
    n = 3.0;
  }
  else if(genus == 4){
    n = 3.92633;
  }
  else if(genus == 5){
    n = 5.25146;
  }
  else{
    n = (double)genus;
  }

  //L = CeilToZZ(sqrt(to_RR(E))*sqrt(n));  
  L = CeilToZZ(to_RR(L)*sqrt(n));

  // Get the basis, {1, rho, omega}.
  init(2*genus);
  basis(G, H);

  GS.G = T.G = B.G = A.G = G; GS.H = T.H = B.H = A.H = H;
  /*
 A.mu0 = ZZ_pX(0, 8) + ZZ_pX(1, 3) + ZZ_pX(2, 7) + ZZ_pX(3, 6) + ZZ_pX(4, 20) + ZZ_pX(5, 26) + ZZ_pX(6, 28) + ZZ_pX(7, 15);
  A.mu1 = ZZ_pX(0, 25) + ZZ_pX(1, 12) + ZZ_pX(2, 4) + ZZ_pX(3, 4) + ZZ_pX(4, 18);
  A.mu2 = ZZ_pX(0, 8) + ZZ_pX(1, 21) + ZZ_pX(2, 19) + ZZ_pX(3, 8) + ZZ_pX(4, 11);
  A.nu0 = ZZ_pX(0, 25) + ZZ_pX(1, 28) + ZZ_pX(2, 9) + ZZ_pX(3, 16) + ZZ_pX(4, 2) + ZZ_pX(5, 15) + ZZ_pX(6, 10) + ZZ_pX(7, 13);
  A.nu1 = ZZ_pX(0, 14) + ZZ_pX(1, 0) + ZZ_pX(2, 8) + ZZ_pX(3, 24) + ZZ_pX(4, 15) + ZZ_pX(5, 25);
  A.nu2 = ZZ_pX(0, 22) + ZZ_pX(1, 1) + ZZ_pX(2, 24) + ZZ_pX(3, 25) + ZZ_pX(4, 13) + ZZ_pX(5, 4);
  A.d   = ZZ_pX(0, 3) + ZZ_pX(1, 24) + ZZ_pX(2, 22) + ZZ_pX(3, 11) + ZZ_pX(4, 17) + ZZ_pX(5, 15) + ZZ_pX(6, 6) + ZZ_pX(7, 7) + ZZ_pX(8,1);
  //A.mu0 = ZZ_pX(0, 9) + ZZ_pX(1, 5) + ZZ_pX(2, 0) + ZZ_pX(3, 2) + ZZ_pX(4, 23) + ZZ_pX(5, 1);
  //A.mu1 = ZZ_pX(0, 11) + ZZ_pX(1, 14) + ZZ_pX(2, 1);// + ZZ_pX(3, 1);
  //A.mu2 = ZZ_pX(0, 17) + ZZ_pX(1, 21) + ZZ_pX(2, 1);// + ZZ_pX(3, 1);
  //A.nu0 = ZZ_pX(0, 26) + ZZ_pX(1, 5) + ZZ_pX(2, 3) + ZZ_pX(3, 11) + ZZ_pX(4, 1);//+ZZ_pX(5, 18);
  //A.nu1 = ZZ_pX(0, 3) + ZZ_pX(1, 1);// + ZZ_pX(2, 26);
  //A.nu2 = ZZ_pX(0, 7) + ZZ_pX(1, 1);// + ZZ_pX(2, 14);
  //A.d   = ZZ_pX(0, 14) + ZZ_pX(1, 1) + ZZ_pX(2, 14) + ZZ_pX(3, 1);// + ZZ_pX(4, 2) + ZZ_pX(5, 1);
  cout<<"A = "<<endl;
  A.print();
  inf_to_cubic(A, T);
  cout<<"1. OK? "<<T.is_valid()<<endl; 
  reduce_ideal(A);
  inf_to_cubic(A, T);
  cout<<"After reduce: 2. OK? "<<T.is_valid()<<endl;
  while(!T.is_valid()){
    T.u++;
  }
  cout<<"T = "<<endl;
  T.print();
  cout<<"3. OK? "<<T.is_valid()<<endl;
  cubic_to_inf(T, A);
  cout<<"A = "<<endl;
  A.print();
  inf_to_cubic(A, T);
  cout<<"T = "<<endl;
  T.print();
  cout<<"4. OK? "<<T.is_valid()<<endl;
  cubic_to_inf(T, A);
  cout<<"nu = "<<A.nu0 + RightShift(A.nu1*rho + A.nu2*omega, prec)<<endl;
 

  A.print();
  
  return;
  */

  reduce_basis(A);
  /*  
  cout<<"Computing near("<<E<<")"<<endl;
  near(E, A);
  baby_step_r1(A, B);
  A.print();
  cout<<"Distance check: "<<A.d0<<" <= "<<E<<" <= "<<B.d0<<endl;

  return;
  */
  long hashsize = NextPrime(4*to_long(L));
  hashbits = NumBits(hashsize);

  ZZVec baby1 = ZZVec(hashsize, q.size()+1);  // First and
  ZZVec baby2 = ZZVec(hashsize, q.size()+1);  // second baby step arrays
  ZZVec baby_deg = ZZVec(hashsize, L.size()+1); // Stores the degrees.
  if((&baby1 == NULL) || (&baby2 == NULL) || (&baby_deg == NULL) ){
    cout<<"Not enough memory available to create baby steps.\n";
    cout<<"Hash table size = "<<hashsize<<".\n"<<"Total memory required = "<<3*hashsize*sizeof(ZZ)<<" bytes.\n\n";
    return;
  }

  //cout<<"Initializing the hashtable."<<endl;

  // Initialize the hashtables.
  for(long i=0; i<hashsize; i++){
    baby_deg[i] = N1;
    clear(baby1[i]);
    clear(baby2[i]);
  }
  //cout<<"Initialized the hashtable."<<endl;
  // Get the basis, {1, rho, omega}.
  init(2*genus);
  basis(G, H);

  cout<<"Making Baby Steps ... \n";
  set = time(NULL);
  order = makeBabySteps1(baby1, baby2, baby_deg, hashsize, GS);
  bstime = time(NULL);
  //for(long i=0; i<hashsize; i++){
  //  cout<<i<<": "<<baby1[i]<<" "<<baby2[i]<<" "<<baby_deg[i]<<endl;
  //}

  if(!IsZero(order)){
    if(order < (E-(genus+2))){
      cout<<"\n"<<"The x-Regulator is R_x = "<<2*order<<".\n";
      cout<<"The S-Regulator is R^S = "<<order<<".\n\n";
    }
    else{
      cout<<"\n"<<"A multiple of R^S is h_0 = "<<order<<".\n"; 
    }

    //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
    cout<<" |h-E|  = "<<abs(order-E)<<endl; 
    //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
    //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
    cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(L32)<<endl;
    cout<<"\n"<<"Baby Steps completed in "<<bstime-set<<" seconds.\n"<<endl;

    return;
  }
  cout<<"\n"<<"Computed "<<L<<" Baby Steps in "<<bstime-set<<" seconds. \nMaking Giant Steps ... \n";

  //cout<<GS.d0<<endl;
  
 
  A = infrastructure_ideal(GS);
  //cout<<"k = "<<0<<" "<<B.d<<" "<<B.is_basis_reduced()<<" "<<B.is_distinguished()<<endl;
  k = 0;
  //cout<<"k = "<<k<<" "<<A.d<<" "<<A.is_basis_reduced()<<" "<<A.is_distinguished()<<endl;

  //return;

  do{
    //A.print();
    //cout<<"k = "<<k<<" "<<A.d<<" "<<A.is_basis_reduced()<<" "<<A.is_distinguished()<<" "<<A.d0<<endl;
    place1 = search(A, baby1, baby2, baby_deg, hashsize);
    
    if(place1 >= 0){
      //cout<<"Found: "<<A.d0<<endl;
      found = 1;
    }
    
    // Search the hash table for the inverse of z.
    inverse(I, A);
    //if(k == 39){
    //  cout<<"k = "<<k<<" "<<I.d<<" "<<I.is_basis_reduced()<<" "<<I.is_distinguished()<<" "<<I.d0<<endl;
    //  I.print();
    //}
    place2 = search(I, baby1, baby2, baby_deg, hashsize);
    
    if(place2 >= 0)
      found = 2;
    
    //if(IsZero(k%percent)){
    //printf("\r");
    //printf("Giant steps are %d %% complete. ", pComplete);
      
    //cout<<flush;
    //pComplete++;
    //}
    
    B = infrastructure_ideal(A);
    A = B*GS;
    s2 = GCD(A.mu2, A.nu2);
    s1 = (A.mu1*A.nu2 - A.nu1*A.mu2)/s2;
    s = A.d;
    //cout<<s<<" "<<s1<<" "<<s2<<endl;
    if(deg(s) != deg(s1))
      baby_step_r1(B, A);
    k++;
      
  }while((place1 < 0) && (place2 < 0));   
  gstime = time(NULL); 
  cout<<"Computed "<<k<<" Giant Steps in "<<gstime-bstime<<" seconds."<<endl;
  //cout<<B.d0<<" "<<baby_deg[place1]<<" place = "<<place1<<endl;
  
  if(found == 1)
    order = (baby_deg[place1] - B.d0)/2;
  else{
    //cout<<"Inverse found it!"<<endl;
    order = (baby_deg[place2] - I.d0)/2;
  }

  cout<<"\n"<<"A multiple of R^S is h_0 = "<<order<<".\n\n"; 
  
  //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<" |h-E|  = "<<abs(order-E)<<endl; 
  //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(L32)<<"\n"<<endl;
  
  // Extract R_S as a factor of h_0.
  R = extract(order, GS.d0/2);
  cout<<endl;
  cout<<"h_0     = "<<order<<endl;
  cout<<"R_S     = "<<R<<endl;
  cout<<"h^*     = "<<order/R<<endl;

  return;
}

// Baby-step, Giant-step algorithm for unit rank 2.

void bsgs2(){
  ZZ order;
  infrastructure_ideal A, A1, A2;
  cubic_ideal B;

  // Get the basis, {1, rho, omega}.
  init(2*genus);
  basis(G, H);

  A1.G = A2.G = A.G = G; A1.H = A2.H = A.H = H; 

  // theta_{1,0}

  A1.mu0 = ZZ_pX(0,3) + ZZ_pX(1, 3) + ZZ_pX(2, 1) + ZZ_pX(3, 4) + ZZ_pX(4, 3) + ZZ_pX(5, 1);
  A1.mu1 = ZZ_pX(0,2) + ZZ_pX(1, 1) + ZZ_pX(2, 6) + ZZ_pX(3, 1); 
  A1.mu2 = ZZ_pX(0,0) + ZZ_pX(1, 2) + ZZ_pX(2, 1);
  A1.nu0 = ZZ_pX(0,1) + ZZ_pX(1, 3) + ZZ_pX(2, 0) + ZZ_pX(3, 4) + ZZ_pX(4, 1);
  A1.nu1 = ZZ_pX(0,2) + ZZ_pX(1, 0) + ZZ_pX(2, 1); 
  A1.nu2 = ZZ_pX(0,3) + ZZ_pX(1, 1);
  A1.d   = ZZ_pX(0,2) + ZZ_pX(1, 5) + ZZ_pX(2, 1);// + ZZ_pX(3, 1);


  A2.mu0 = ZZ_pX(0,3) + ZZ_pX(1, 6) + ZZ_pX(2, 3) + ZZ_pX(3, 4) + ZZ_pX(4, 2) + ZZ_pX(5, 1);
  A2.mu1 = ZZ_pX(0,1) + ZZ_pX(1, 5) + ZZ_pX(2, 5) + ZZ_pX(3, 1); 
  A2.mu2 = ZZ_pX(0,1) + ZZ_pX(1, 1) + ZZ_pX(2, 1);
  A2.nu0 = ZZ_pX(0,4) + ZZ_pX(1, 4) + ZZ_pX(2, 2) + ZZ_pX(3, 5) + ZZ_pX(4, 1);
  A2.nu1 = ZZ_pX(0,6) + ZZ_pX(1, 3) + ZZ_pX(2, 1); 
  A2.nu2 = ZZ_pX(0,2) + ZZ_pX(1, 1);
  A2.d   = ZZ_pX(0,1) + ZZ_pX(1, 3) + ZZ_pX(2, 4) + ZZ_pX(3, 1);


  
  A = A1*A2;

  A.print();
  //reduce_basis(A);

  //A.print();

  //baby_step_0_r2(A, A1);
  //A1.print();
  
  inf_to_cubic(A, B);
  cout<<"Have not implemented Baby-step, Giant-step for unit rank 2.\n"<<"In fact, nobody has, anywhere.\n\n";
  return;
}


// Returns 1 if a is a cube in ZZ_p.
// Returns 0 if not.

int isCube(ZZ_p a){
  if(q%3 == 2)
    return(1);
  else
    return(IsOne(power(a, (q-1)/3))); 
}

// This function reads in the data from the file input.

int getinput(char *input){
  FILE *inputptr;
  char thisLine[128];
  char token[64], value[64];
  int dD, dG, dH, rand;
  int coefg[10];
  int coefh[10];
  int i;
  vec_pair_ZZ_pX_long factors;

  maxCollision = 0;

  SetSeed( to_ZZ( time(NULL) ) );

  inputptr = fopen(input, "r");

  if(inputptr == NULL){
    cout<<"Error: Could not open "<<input<<".\n\n";
    return(0);
  }

  while(!feof(inputptr)){
    thisLine[0] = 0;
    fgets(thisLine, 128, inputptr);
    // Skip comments and blank space.
    while((thisLine[0] =='#') || isspace(thisLine[0])){
      fgets(thisLine, 128, inputptr);
    }
    sscanf(thisLine, "%s%s", token, value);
    if(strncmp(token, "q:", 2)==0){
      q = to_ZZ(value);
      // Set the characterisitc.
      if(q<=to_ZZ(3)){
	cout<<"Error: Characteristic q must be greater than 3.\n\n";
	return(0);
      }
      ZZ_p::init(to_ZZ(q));
    }
    else if(strncmp(token, "l:", 2)==0)
      lambda = atoi(value);
    else if(strncmp(token, "bsgs:", 5)==0)
      bsgs = atoi(value);
    else if(strncmp(token, "kang:", 5)==0)
      kang = atoi(value);
    else if(strncmp(token, "degreeG:", 8)==0)
      dG = atoi(value);
    else if(strncmp(token, "degreeH:", 8)==0)
      dH = atoi(value);
    else if(strncmp(token, "random:", 7)==0){
      rand = atoi(value);
      if( rand == 2 )
	irred = 1;
      else
	irred = 0;
    }
    else if((strncmp(token, "g0:", 3)==0) && !rand) {
      coefg[0] = atoi(value);
      for(i=1; i<=dG; i++){
	fgets(thisLine, 128, inputptr);
	sscanf(thisLine, "%s%s", token, value);
	coefg[i] = atoi(value);
      }

      // Make the polynomial G(x).
      for(i=0; i<=dG;i++)
	G += ZZ_pX(i,coefg[i]);
    }
    else if((strncmp(token, "h0:", 3)==0) && !rand) {
      // We set H(x) = x if deg(H) = 1 since that simplifies the arithmetic.
      if(dH == 0){
	H = ZZ_pX(0,1);
	D = G;
      }
      else if(dH == 1){
	SetX(H);
	D = LeftShift(G, 2);
      }
      else{
	coefh[0] = atoi(value);
	for(i=1; i<=dH; i++){
	  fgets(thisLine, 128, inputptr);
	  sscanf(thisLine, "%s%s", token, value);
	  coefh[i] = atoi(value);
	}
	// Make the polynomial H(x).
	for(i=0; i<=dH;i++)
	  H += ZZ_pX(i,coefh[i]);

	D = G*sqr(H);
      }
    }
    else{}
  }

  // Generate random G and H if desired.
  if(rand) {
    G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
    if(irred){
      while(!DetIrredTest(G))
	G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
    }

    // Make sure G is squarefree.
    while(!IsOne(GCD(G, diff(G)))){
      G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
    }
    if(dH == 0){
      H = ZZ_pX(0,1);
      D = G;
    }
    else if(dH == 1) {
      SetX(H);
      
      // Make sure that G and H are relatively prime.
      // And that G is squarefree.
      if(IsZero(coeff(G,0)) || !IsOne(GCD(G, diff(G))) ){
	while(IsZero(coeff(G,0)) || !IsOne(GCD(G, diff(G))) )
	  G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
      }
      D = LeftShift(G, 2);
    }
    else {
      H = ZZ_pX(dH, 1) + random_ZZ_pX(dH);
      if(irred){
	do{
	  while(!DetIrredTest(H))
	    H = ZZ_pX(dH, 1) + random_ZZ_pX(dH);
	}while(G == H);
      }

      // Make sure G and H are relatively prime and G is squarefree.
      while(!IsOne(GCD(G,H)) || !IsOne(GCD(G, diff(G)))  ){
	H = ZZ_pX(dH, 1) + random_ZZ_pX(dH);
      }

      D = G*sqr(H);
    }
    
  }
  else{
    // Make sure that (G,H) = 1 and G is squarefree.
    if( !IsOne(GCD(G,H)) || !IsOne(GCD(G, diff(G))) ){
      cout<<"Generating new G(x). Original was either not squarefree or not coprime to H(x).\n\n";
    }
    while( !IsOne(GCD(G,H))  || !IsOne(GCD(G, diff(G))) ){
      G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
    }
    D = G*sqr(H);
  }

  // It may be that G and H are irreducible anyway.
  if(!irred){
    if(dH >=2){
      if(DetIrredTest(H) && DetIrredTest(G))
	irred = 1;
    }
    else{
      if(DetIrredTest(G))
	irred = 1;
    }
  }

  // Find the numbers and degrees of the factors of G and H.
  for(i=0; i<9; i++)
    factordegs[i]=0;
  if(irred){
    factordegs[dG-1]++;
    if(dH)
      factordegs[dH-1]++;
  }
  else{
    if(dH==1)
      factordegs[0]++;
    else{
      factors = berlekamp(inv(LeadCoeff(H))*H);
      for(i=0; i<factors.length(); i++){
	factordegs[deg(factors[i].a)-1]++;
      }    
    }
    factors = berlekamp(inv(LeadCoeff(G))*G);
    
    for(i=0; i<factors.length(); i++){
      //cout<<factors[i].a<<endl;
      factordegs[deg(factors[i].a)-1]++;
    } 
  }
  // Set the genus and unit rank.
  if ((dG + 2*dH)%3) {
    rank = 0;
    genus = dG + dH - 1;
  }
  else {
    genus = dG + dH - 2;
    if(q%3 == 1)
      rank = 2;
    else if(isCube(LeadCoeff(D)))
      rank = 1;
    else
      rank = 0;
  }

  if(genus == 3){
    alpha = 0.27183722; //0.26864405;
  }
  else if(genus == 4){
    alpha = 0.19169713; //0.19208172; //0.19336089; //0.17759;
  }
  else if(genus == 5){
    alpha = 0.17981087;
  }
  else if(genus == 6){
    alpha = 0.15978348; //0.15841058;
  }
  else if(genus == 7){
    alpha = 0.10909269;
  }
  else{
    // We have no experimental data for this case, so
    // we're winging it here.
    alpha = 1.0/(2.0*(genus-1.0));
  }

  // In the irreducible case, h = 1 (mod 3).
  // Otherwise h = 0, 1 (mod 3)
  // Taking advantage of this will lessen the number of steps.
  if(rank == 0){
    if(irred && (dH == 0)){
      alpha/=3.0;
    }
    else{
      alpha*=(2.0/3.0);
    }
  }
  fclose(inputptr);
 
  // For hash storage. 

  if(genus < 5)
    type = 1;
  else
    type = 2;

  return(1);
 
}

void printInitData(){
  int dD = deg(D), dG = deg(G), dH = deg(H);
  int i;

  cout<<"\n";
  cout<<"Working with the function field F_"<<q<<"(x, y), where"<<endl;
  cout<<"y^3 = G(x) H^2(x) = ";

  if (dD > 1){
    if(!IsOne(coeff(D,dD)))
      cout<<coeff(D,i);
    cout<<"x^"<<dD;
  }
  for(i=dD-1; i>1; i--){
    if(!IsZero(coeff(D, i))) {
      cout<<" + ";
      if(!IsOne(coeff(D,i)))
	cout<<coeff(D,i);
      cout<<"x^"<<i;
    }
  }
  if(!IsZero(coeff(D,1))){
    if(dD > 1) { cout<<" + "; }
    if(!IsOne(coeff(D,1)))
      cout<<coeff(D,1);
    cout<<"x";
  }
  if(!IsZero(coeff(D,0)))
    cout<<" + "<<coeff(D,0)<<endl;
  else
    cout<<"\n";

  cout<<"G(x) = ";
  if(dG > 1){
    if(!IsOne(coeff(G,dG)))
      cout<<coeff(G,i);
    cout<<"x^"<<dG;
  }
  for(i=dG-1; i>1; i--){
    if(!IsZero(coeff(G, i))) {
      cout<<" + ";
      if(!IsOne(coeff(G,i)))
	cout<<coeff(G,i);
      cout<<"x^"<<i;
    }
  }
  if(!IsZero(coeff(G,1))){
    if(dG > 1) {cout<<" + ";}
    if(!IsOne(coeff(G,1)))
      cout<<coeff(G,1);
    cout<<"x";
  }
  if(!IsZero(coeff(G,0)))
    cout<<" + "<<coeff(G,0)<<endl;
  else
    cout<<"\n";

  cout<<"H(x) = ";
  if(dH == 0){
    cout<<"1\n";
  }
  else  if(dH == 1){
    cout<<"x\n";
  }
  else{
    if(!IsOne(coeff(H,dH)))
      cout<<coeff(H,i);
    cout<<"x^"<<dH;
  
    for(i=dH-1; i>1; i--){
      if(!IsZero(coeff(H, i))) {
	cout<<" + ";
	if(!IsOne(coeff(H,i)))
	  cout<<coeff(H,i);
	cout<<"x^"<<i;
      }
    }
    if(!IsZero(coeff(H,1))){
      if(dH > 1){ cout<<" + "; }
      if(!IsOne(coeff(H,1)))
	cout<<coeff(H,1);
      cout<<"x";
    }
    if(!IsZero(coeff(H,0)))
      cout<<" + "<<coeff(H,0)<<endl;
  }
  
  cout<<"The curve has UNIT RANK = "<<rank<<endl;
  cout<<"                  GENUS = "<<genus<<"\n\n";

}
void printInitDataFile(char *file){
  FILE *outputptr;
  int dD = deg(D), dG = deg(G), dH = deg(H);
  int i;

  // For now, all this prints out to the screen.
  cout<<"Working with the function field F_"<<q<<"(x, y), where"<<endl;
  cout<<"y^3 = G(x) H^2(x) = ";

  if (dD > 1){
    if(!IsOne(coeff(D,dD)))
      cout<<coeff(D,i);
    cout<<"x^"<<dD;
  }
  for(i=dD-1; i>1; i--){
    if(!IsZero(coeff(D, i))) {
      cout<<" + ";
      if(!IsOne(coeff(D,i)))
	cout<<coeff(D,i);
      cout<<"x^"<<i;
    }
  }
  if(!IsZero(coeff(D,1))){
    if(dD > 1) { cout<<" + "; }
    if(!IsOne(coeff(D,1)))
      cout<<coeff(D,1);
    cout<<"x";
  }
  if(!IsZero(coeff(D,0)))
    cout<<" + "<<coeff(D,0)<<endl;
  else
    cout<<"\n";

  cout<<"G(x) = ";
  if(dG > 1){
    if(!IsOne(coeff(G,dG)))
      cout<<coeff(G,i);
    cout<<"x^"<<dG;
  }
  for(i=dG-1; i>1; i--){
    if(!IsZero(coeff(G, i))) {
      cout<<" + ";
      if(!IsOne(coeff(G,i)))
	cout<<coeff(G,i);
      cout<<"x^"<<i;
    }
  }
  if(!IsZero(coeff(G,1))){
    if(dG > 1) {cout<<" + ";}
    if(!IsOne(coeff(G,1)))
      cout<<coeff(G,1);
    cout<<"x";
  }
  if(!IsZero(coeff(G,0)))
    cout<<" + "<<coeff(G,0)<<endl;
  else
    cout<<"\n";

  cout<<"H(x) = ";
  if(dH > 1){
    if(!IsOne(coeff(H,dH)))
      cout<<coeff(H,i);
    cout<<"x^"<<dH;
  }
  for(i=dH-1; i>1; i--){
    if(!IsZero(coeff(H, i))) {
      cout<<" + ";
      if(!IsOne(coeff(H,i)))
	cout<<coeff(H,i);
      cout<<"x^"<<i;
    }
  }
  if(!IsZero(coeff(H,1))){
    if(dH > 1){ cout<<" + "; }
    if(!IsOne(coeff(H,1)))
      cout<<coeff(H,1);
    cout<<"x";
  }
  if(!IsZero(coeff(H,0)))
    cout<<" + "<<coeff(H,0)<<endl;
  else
    cout<<"\n";
  
  cout<<"The curve has UNIT RANK = "<<rank<<endl;
  cout<<"                  GENUS = "<<genus<<"\n\n";
  
}

// Return 1 if the order is probably correct.
// Return 0 if it is not.

int check(ZZ &order){
  cubic_ideal g, temp;

  if(checkPrimes(order)){
    g.G = G; g.H = H;
    temp.G = G; temp.H = H;
    
    for(int i = 0; i < 20; i++){
      random(g);
      temp = g^order;
      // This isn't the order. Keep searching.
      if(!IsOne(temp.s)){
	cout<<"\nFound the wrong order.\n\n";
	return 0;
      }
    }

    return 1;
  }
  cout<<"\nFound the wrong order.\n\n";
  return 0;
}

int main(int argc, char *argv[]){
  ZZ factor;
  time_t set, time1, time2, time3;

  if(!((argc == 2) || (argc == 3))){
    cout<<"Usage: classnumber inputfile [outputfile]\n\n";
    return 0;
  }
  else {
    if(!getinput(argv[1])) {
      return(0);
    }
  }

  // Print the initial data.
  if(argc == 2)
    printInitData();
  else
    printInitDataFile(argv[2]);
  /*
  cubic_ideal x,y,z;
  ZZ e = to_ZZ(14274);

  x.G = y.G = z.G = G;
  x.H = y.H = z.H = H;
    
  x.s = ZZ_pX(6,1)+ZZ_pX(5,0)+ZZ_pX(4,0)+ZZ_pX(3,5)+ZZ_pX(2,7)+ZZ_pX(1,8)+ZZ_pX(0,7);
  y.s2 = x.s2 = x.s1 = ZZ_pX(1,0)+ZZ_pX(0,1);
  x.u = ZZ_pX(5,2)+ZZ_pX(4,7)+ZZ_pX(3,0)+ZZ_pX(2,8)+ZZ_pX(1,6)+ZZ_pX(0,5);
  x.v = ZZ_pX(5,5)+ZZ_pX(4,1)+ZZ_pX(3,0)+ZZ_pX(2,8)+ZZ_pX(1,2)+ZZ_pX(0,11);
  x.w = ZZ_pX();

  y.s = ZZ_pX(5,1)+ZZ_pX(4,4)+ZZ_pX(3,6)+ZZ_pX(2,5)+ZZ_pX(1,2)+ZZ_pX(0,2);
  y.s1 = ZZ_pX(1,1)+ZZ_pX(0,6);
  y.u = ZZ_pX(3,8)+ZZ_pX(2,7)+ZZ_pX(1,8)+ZZ_pX(0,5);
  y.v = ZZ_pX(4,8)+ZZ_pX(3,4)+ZZ_pX(2,6)+ZZ_pX(1,6)+ZZ_pX(0,0);
  y.w = ZZ_pX(0,8);

  

  cout<<"Valid ideals? "<<x.is_valid()<<" "<<y.is_valid()<<endl;
  z = x*y;
  //z = x^e;
  //square(z,x);
  cout<<"Valid result? "<<z.is_valid()<<endl;

  z.print();
  return 0;
  */  
  // Step 1: Compute an approximation E of the class number h
  //         and an integer L such that |h - E| < L^2.
  set = time(NULL);

  approxh();

  time1 = time(NULL);

  cout<<"Phase 1: Completed in "<<time1-set<<" seconds.\n\n";

  // Step 2: Use extra information about h in (E-L^2, E+L^2),
  //         such as its distribution in the interval
  //         or h mod r for small primes r.

  // Find elements of small order.
  if(rank == 0){
    factor = smallOrder0();
    if(IsOne(factor))
      cout<<"Phase 2: h = 1 (mod 3).\n\n";
    else
      cout<<"Phase 2: h = 0 (mod "<<factor<<").\n\n";
    //cout<<"Phase 2: Completed in "<<time2-time1<<" seconds.\n\n";
  }

  time2 = time(NULL);

  // Step 3: Find h in the interval (E-L^2, E+L^2) via
  //         A) Baby Step, Giant Step or
  //         B) Pollard's Kangaroo
  if(rank == 0){
    if(bsgs)
      bsgs0(factor);
    if(kang)
      kangaroo0(factor);
  }
  else if(rank == 1)
    bsgs1();
  else
    bsgs2();
  time3 = time(NULL);
  cout<<"Phase 3: Completed in "<<time3-time2<<" seconds.\n\n";
  cout<<"Total running time: "<<time3-set<<" seconds.\n\n";

}
