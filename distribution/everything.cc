#include "everything.h"

NTL_CLIENT


// Finds an approximation E of h and a bound L such that
// |h-E|<L^2.

void approxh(){
  int s1, s2;
  int n, v, divlim; // v|n
  int i, j;
  RR B;  // Set logE2 = A(K) + B.
  ZZ B1; // Used to find the trivial bound E1.
  RR C, S; // Temp variables.
  RR Q = to_RR(q);
  RR A;
  ZZ SvSum;

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
 
  /***********************************************/
  /* Precompute the S_v(j) for 1 <= v <= lambda. */
  /***********************************************/

  // This will most likely be the longest running loop 
  // before the baby-step, giant-step portion.
  for(v=1; v<=lambda; v++){
    if((q%3 == 1)||(v%2 == 0)) {
      Sv1(v, 1);
      Sv1(v, 3);
    } else{

      // Note: This must be fixed for q=2(mod 3)
      // since S_v(a) takes on four values in this case:
      // S_v(1), S_v(2), S_v(3), and S_v(6).
      Sv2(v, 1);
      Sv2(v, 2);
    }
  }

  //cout<<SvCache[0]<<" "<<SvCache[1]<<" "<<SvCache[2]<<" "<<SvCache[3]<<endl;

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

  /************************************************/
  /* Approximate the square root of the error, L. */
  /************************************************/

  if(q%3 == 1){

    /*********************************************/
    /* These are the trivial bounds: E1 and L1^2 */
    /*********************************************/
    E1r = exp(logE1);
    Es[0] = RoundToZZ(E1r);
    L2[0] = CeilToZZ(E1r*expm1(psi1)+1.0/2.0);
    Ls[0] = CeilToZZ(sqrt(alpha1*to_RR(L2[0])));

    /*************************************/
    /* These are the bounds: E2 and L2^2 */
    /*************************************/
    	
    // We're taking advantage of easy inversions here.
    E2r = exp(logE2);
    Es[1] = RoundToZZ(E2r);
    L2[1] = CeilToZZ(E2r*expm1(psi2)+1.0/2.0);
    Ls[1] = CeilToZZ(sqrt(alpha2*to_RR(L2[1])));

    /*************************************/
    /* These are the bounds: E3 and L3^2 */
    /*************************************/

    // E3 = E2.

    // Finish assembling psi_3.

    SvSum = ZZ();
    i = lambda + 1;
    // Finding |SUM_{v|(lambda+1), v != lambda +1}vS_v((lambda+1)/v)|
    for(j=1; j<i; j++){
      if(j%i == 0){
	if((i/j)%3)
	  SvSum += j*SvCache[2*(j-1)];
	else
	  SvSum += j*SvCache[2*j-1];
      }
    }

    psi3 += to_RR(abs(SvSum))/(i*power(Q, i));
    Es[2] = Es[1];
    L2[2] = CeilToZZ(E2r*expm1(psi3)+1.0/2.0);
    Ls[2] = CeilToZZ(sqrt(alpha3*to_RR(L2[2])));

    //RoundToZZ(E, E2);
    //L1 = CeilToZZ(sqrt(alpha*to_RR(L12)));
    //L2 = CeilToZZ(sqrt(alpha*to_RR(L22)));
    //L3 = CeilToZZ(sqrt(alpha*to_RR(L32)));

    // Find the best interval to use.
    /*
      if(E1-L12 > E2-L32){ 
      if(E1+L12 < E2+L32)
      best = 1;
      else{
      if(L12 < L32){
      E = E1;
      L = CeilToZZ(sqrt(alpha*to_RR(L12)));
      best = 0;
      }
      else{
      E = E2;
      L = CeilToZZ(sqrt(alpha*to_RR(L32)));
      best = 0;
      }
      }
      }
      else{
      if(E1+L12 > E2+L32)
      best = 3;
      else{
      if(L12 < L32){
      E = E1;
      L = CeilToZZ(sqrt(alpha*to_RR(L12)));
      best = 0;
      }
      else{
      E = E2;
      L = CeilToZZ(sqrt(alpha*to_RR(L32)));
      best = 0;
      }
      }
      }
      if(best == 1){
      E = E1;
      L = CeilToZZ(sqrt(alpha*to_RR(L12)));
      }
      if(best == 3){
      E = E2;
      L = CeilToZZ(sqrt(alpha*to_RR(L32)));
      }
    */
  }
  else{
    if(deg(D)%3 != 0) j=0;
    else j=2;
    C = power(Q,lambda);
    for(i=lambda + 1; i<= 2*lambda; i++){
      C*=Q;
      if(i%2 == 0)
	B+=(j + 2.0*genus*power(Q,i/2) + (2.0*Q/(Q-1))*(power(Q,i/2) - 1))/(i*C);
    }
  }
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
	  L+=2;
	}
	else{
	  // This is the inert case.
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
	  L+=2;
	}
	else{
	  // This is the inert case.
	  L--;
	}
	j++;
      }
    }
  }
  // We use a nice trick to run through every irreducible 
  // degree 2 polynomial.
  else if (v == 2){
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
	      L+=2;
	    }
	    else{
	      // This is the inert case.
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
	      L+=2;
	    }
	    else{
	      // This is the inert case.
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
    if(factordegs[v-1]){
      while(IsOne(LeadCoeff(P))){
	// Formerly if(ProbIrredTest(P, iter=1)){
	if(DetIrredTest(P)){
	  if(IsZero(D%P));
	  else if(chi(P) == 1){
	    L+=2;
	  }
	  else{
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
	    L+=2;
	  }
	  else{
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

  L*=2;
  SvCache[loc] = L;
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

  GCD(T, F, P);
  if(!IsOne(T))
    return ZZ_p();

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
  ZZ k=to_ZZ(3), factor;
  vec_pair_ZZ_pX_long factors;

  factors = berlekamp(inv(LeadCoeff(G))*G);
  if(IsOne(H))
    power(factor, k, factors.length()-1);
  else
    power(factor, k, factors.length());

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


// Baby-step, Giant-step algorithm for unit rank 0 with the first approximation.
// Input: factor - an integer that we know is a factor of the order h.
// 
// Knowing so speeds up the algorithm by a factor of sqrt(so).

void bsgs0(ZZ& factor, int run){
  int found = 0, genTrials=0, test1=200, test2=20, i;
  long hashsize, place1, place2;
  cubic_ideal g; // A subgroup generator.
  cubic_ideal z, temp, zinv;
  cubic_ideal gs; // The giant step.
  cubic_ideal I; // The identity.
  ZZ start, end=E+L2[run];
  ZZ order, order2, giantstep;
  ZZ ordg;
  ZZ numSteps; 
  ZZ stepLen = factor;
  ZZ k = ZZ(), r;
  double n;   // n = #time for one giant step/#time for one baby step.

  if( genus == 3 ){
    n = 1.38;
  }
  else if( genus == 4 ){
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

  order0 = ZZ();
  order1 = ZZ();

  z.G = G;
  z.H = H;
  temp.G = G;
  temp.H = H;
  gs.G = G;
  gs.H = H;

  E = Es[run];
  L = Ls[run];

  // The number of Baby Steps. It will be much nicer if they are even.
  numSteps = CeilToZZ(to_RR(L)*sqrt(n)/sqrt(to_RR(stepLen)));
  numSteps += (numSteps%2);
  
  hashsize = NextPrime(4*to_long(numSteps));
  hashbits = NumBits(hashsize);

  // We're only going to store a couple coefficients of s 
  // from the ideal to save time and space.
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

  do {
    for(long i=0; i<hashsize; i++){
      baby_loc[i] = N1;
      //baby[i] = I;
      clear(baby1[i]);
      clear(baby2[i]);
    }

    // Get a random (sub)group generator.
    random(g);

    //g.s = ZZ_pX(0,16)+ZZ_pX(1, 6)+ZZ_pX(2, 1)+ZZ_pX(3, 0)+ZZ_pX(4, 0)+ZZ_pX(5, 0);
    //g.s1 = gs.s2 = ZZ_pX(0,1);
    //g.u = ZZ_pX(0,10)+ZZ_pX(1, 17)+ZZ_pX(2, 0)+ZZ_pX(3, 0)+ZZ_pX(4, 0);
    //g.v = ZZ_pX(0,11)+ZZ_pX(1, 7)+ZZ_pX(2, 0)+ZZ_pX(3, 0)+ZZ_pX(4, 0);
    //g.w = ZZ_pX();

    //cout<<"Generator g = "<<endl;
    //g.print();
    
    setb[run] = time(NULL);
    ordg = makeBabySteps(g, baby1, baby2, baby_loc, hashsize, numSteps, stepLen, 1);

    //cout<<"Baby steps completed."<<endl;
    
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
	h = ordg;
	return;
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
  timeb[run] = time(NULL);
  genTrials = 0;
  babysteps[run] = numSteps;

  /***************/
  /* GIANT STEPS */
  /***************/

  //cout<<"Beginning giant steps."<<endl;

  // Set the giant step distance and the first giant step.

  if(irred && (degH == 0)){
    giantstep = stepLen*numSteps;
    gs = g^giantstep;
  }
  else{
    giantstep = stepLen*numSteps/2;
  }
  gs = g^giantstep;
  z = cubic_ideal(gs);

  //gs.s = ZZ_pX(0,1)+ZZ_pX(1, 0)+ZZ_pX(2, 1)+ZZ_pX(3, 7)+ZZ_pX(4, 15)+ZZ_pX(5, 1);
  //gs.s1 = gs.s2 = ZZ_pX(0,1);
  //gs.u = ZZ_pX(0,18)+ZZ_pX(1, 10)+ZZ_pX(2, 5)+ZZ_pX(3, 2)+ZZ_pX(4, 6);
  //gs.v = ZZ_pX(0,10)+ZZ_pX(1, 3)+ZZ_pX(2, 13)+ZZ_pX(3, 3)+ZZ_pX(4, 11);
  //gs.w = ZZ_pX();

  //cout<<"generator:"<<endl;
  //g.print();
  //z = cubic_ideal(gs);

  //cout<<"Giant step = "<<endl;
  //z.print();

  setg[run] = time(NULL);

  do{
    do{
      //cout<<"Searching "<<endl;
      place1 = search(z, baby1, baby2, baby_loc, hashsize);      
    
      //cout<<"place1 = "<<place1<<endl;

      if(place1 >= 0)
	found = 1;
      
      // Search the hash table for the inverse of z.

      //cout<<"z inverse: "<<endl;
      inverse(temp, z);
      //cout<<"temp"<<endl;
      reduce(zinv, temp, D);
      //cout<<"reduced:"<<endl;
      //zinv.print();
      
      place2 = search(zinv, baby1, baby2, baby_loc, hashsize);
      
      if(place2 >= 0)
	found = 2;

      temp = cubic_ideal(z);

      z = temp*gs;
      k++;
    }while((place1 < 0)&&(place2<0));
    
    // At this point we've found a match.
    start = E - (E%stepLen);
    if(place1 >= 0){
      order = start - giantstep*k + baby_loc[place1];
    }
    else{
      order = start + giantstep*k + baby_loc[place2];
    }

    timeg[run] = time(NULL);

    // Check to see if it is valid.
    if(checkPrimes(order)){

      // Test several elements to see if this is the order.
      genTrials = 0;
      while(genTrials < test2){
	random(g);	
	temp = g^order;
	// This isn't the order. Keep searching.
	if(temp != I){
	  found = 0;
	  break;
	}
	else
	  genTrials++;
      }
    }
    else
      found = 0;

    // If two matches were found, then check the other match.
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
	    found = 0;
	    break;
	  }
	  else
	    genTrials++;
	}
      }
      else
	found = 0;
    }
  } while(!found);
  
  giantsteps[run] = k;

  h = order;

  baby1.kill();
  baby2.kill();
  baby_loc.kill();
 
  return;
}

// Input: g: A (sub)group generator.
//        baby[]: A hash table to store baby steps.
//        loc[]: An array that gives the power of g stored in the array.
//        hashsize: the size of the hashtable.
//        steps: number of baby steps to compute.
//        stepLen: The length of a baby step.
// Return: either the order of g, or 0 if the identity was not encountered. 

ZZ makeBabySteps(cubic_ideal &g, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize, ZZ& steps, ZZ& stepLen, int verbose){
  // The identity and a random (sub)group generator.
  cubic_ideal h, z, step;
  long j, index;
  ZZ i = to_ZZ(1), start, steps0;
  
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
  h = g^start;
  step = g^stepLen;

  // If D(x) is irreducible, h = 1 (mod 3). 
  // Work in this equivalence class.
  if(irred && (degH == 0)){
    start++;

    z = h;
    h = z*g;
    
    // Work in the equivalence classes 1 mod stepLen.
    
    for(i=0; i<steps; i++){
      // If the baby step is the identity, then start + i*stepLen + 1
      // is a multiple of the order of g, but not necessarily the
      // group order.
      // Check small primes, especially those = 2 mod 3.
      if(IsOne(h.s)){
	if(!IsZero(order0))
	  order1 = order0;
	order0 = start+i*stepLen + 1;
	if(checkPrimes(order0)){
	  return (order0);
	}
      }
      // Otherwise, we insert the value into the hash table.
      insert(h, baby1, baby2, loc, stepLen*i+1, hashsize); 

      // Make the next baby step.
      z = cubic_ideal(h);
      h = z*step;
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
    
    for(i=0; i<steps0; i++){
    
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
      
      // Otherwise, we insert the value into the hash table.
      insert(h, baby1, baby2, loc, stepLen*i, hashsize); 

      // Make the next baby step.
      z = cubic_ideal(h);
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
      h1 = z*step;
    }
  }

  return ZZ();
}

// The hash function.

long hash(cubic_ideal &A, long size){
  long value;
  int i, d = deg(A.u);
  
  if(IsZero(A.u)){
    if(IsZero(A.v)){
      if(IsOne(A.s))
	return 1;
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
}
// Insert an ideal into the baby step hash table.

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

// Returns the position of A in the hash table, or -1 if it is not found.

inline long search(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize){
  long j, index;
  ZZ bvalue1, bvalue2;
  int i;

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
    index+=(j*j);
    index%=hashsize;
    j++;
    if(loc[index] == N1)
      return(-1);
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


// Baby-step, Giant-step algorithm for unit rank 1.

//void bsgs1(){
//  cout<<"Have not implemented Baby-step, Giant-step for unit rank 1.\n\n";
//  return;
//}

// Baby-step, Giant-step algorithm for unit rank 2.

//void bsgs2(){
//  cout<<"Have not implemented Baby-step, Giant-step for unit rank 2.\n"<<"In fact, nobody has, anywhere.\n\n";
//  return;
//}


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
  int dD;
  int i;

  maxCollision = 0;

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
    else if(strncmp(token, "degreeG:", 8)==0)
      degG = atoi(value);
    else if(strncmp(token, "degreeH:", 8)==0)
      degH = atoi(value);
    else if(strncmp(token, "irreducible:", 12)==0)
      irred = atoi(value);
    else if(strncmp(token, "samples:", 8)==0)
      samples = atoi(value);
    else{}
  }

  fclose(inputptr);
  return(1);
}

// Set other global constants that will be the same for 
// each polynomial: rank, genus, lambda.

void setConstants(){
  int i, j, l, n;
  RR Q = to_RR(q);  
  RR C, S, T;

  /********************************/
  /* Set the genus and unit rank. */
  /********************************/

  if ((degG+2*degH)%3) {
    rank = 0;
    genus = degG + degH -1;
  }
  else {
    genus = degG + degH - 2;
    if(q%3 == 1)
      rank = 2;
    //else if(isCube(LeadCoeff(D)))
    //  rank = 1;
    else
      rank = 0;
  }

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

  /**************/
  /* Set psi_1. */
  /**************/

  T = sqrt(Q);

  psi1 = log(T) - log(T-1.0);
  for(i=1; i<=lambda; i++){
    psi1 -= inv(i*power(T,i));
  }
  psi1*=(2*genus);
  psi1 += 2*(log(Q) - log(Q-1));
  for(i=1; i<=lambda; i++){
    psi1 -= (2*inv(i*power(Q,i)));
  }

  /**************/
  /* Set psi_2. */
  /**************/

  i = lambda+1;
  n = i+1;
  
  psi2 = 2.0*genus/(i*power(T,i)) + ((2*genus + 4)*T*Q*power(T, -n))/(n*(T-1.0)*(Q-1.0)) + 2.0/(i*power(Q, i));
  
  if(i%2)
    psi2 += (2.0*(Q/(Q-1.0))*power(Q, -i)*(power(Q,i/3)-1.0)/i);
  else
    psi2 += (2.0*(Q/(Q-1.0))*power(Q, -i)*(power(Q,i/2)-1.0)/i);

  /********************************************/
  /* Set psi_3, except for the vS_v(n/v) sum. */
  /********************************************/

  // This is a little bit of a speed-up.
  // Instead of dividing odd i by 3 and even i by 2,
  // we divide i by its smallest prime factor.
  // There isn't a speed-up until i=5 is reached,
  // i.e. unless lambda > 2.
    
  for(j=0; j<15; j++){
    if(i%primes[j] == 0){
      l = primes[j];
      break;
    }
  }
  S = 1.0/to_RR(l);
  C = pow(Q, S);
  S = Q/C;
  
  // Finding |SUM(x_j^{lambda+1}|
  if(deg(D)%3)
    psi3 = 2.0/(i*power(Q, i));
  else
    psi3 = RR();
  
  psi3 += (2.0*genus/(i*power(T,i)) + (2*genus*T*Q*power(T, -n))/(n*(T-1)*(Q-1)) + 4.0*(Q/(Q-1))*(S/(S-1))*power(S, -n)/n);
  

  
  /**********************************/
  /* Set the Hasse Interval bounds. */
  /**********************************/

  HasseLow = CeilToZZ(power(T - 1, 2*genus));
  HasseHigh= TruncToZZ(power(T + 1, 2*genus));

  /*************/
  /* Set alpha */
  /*************/

  if(genus == 3){        // For q = 10009.
    alpha1 = 0.27223307; //0.27031818;
    alpha2 = 0.20405624; //0.20234914;
    alpha3 = 0.27183722; //0.26906175;
  }
  else if(genus == 4){
    alpha1 = 0.19236281; 
    alpha2 = 0.15365800; //0.15560498;
    alpha3 = 0.19169713; //0.19336089;
  }
  else if(genus == 5){ // From q = 97.
    alpha1 = 0.18195632;
    alpha2 = 0.17143328;
    alpha3 = 0.17981087;
  }
  else if(genus == 6){
    alpha1 = 0.15997123;
    alpha2 = 0.15679484;
    alpha3 = 0.15978348;
  }
  else if(genus == 7){
    alpha1 = 0.11428348;
    alpha2 = 0.10135344;
    alpha3 = 0.10909269;
  }
  else{
    alpha2 = 1.0/(2.0*genus-1.0);
    alpha1 = alpha3 = alpha2*4.0/3.0;
  }
  
  // In the irreducible case, h = 1 (mod 3).
  // Otherwise h = 0, 1 (mod 3)
  // Taking advantage of this will lessen the number of steps.
  if(irred && (degH == 0)){
    alpha1/=3.0;
    alpha2/=3.0;
    alpha3/=3.0;
  }
  else{
    alpha1*=(2.0/3.0);
    alpha2*=(2.0/3.0);
    alpha3*=(2.0/3.0);
  }
}

// Get a random polynomial for testing.
// G and H are both monic. If deg(H) = 1, then H(x) = x.

void getPolynomial(){
  int i, degD;
  vec_pair_ZZ_pX_long factors;

  // Get a random degree degG polynomial.

  G = ZZ_pX(degG, 1) + random_ZZ_pX(degG);

  if(irred){
    while(!DetIrredTest(G))
      G = ZZ_pX(degG, 1) + random_ZZ_pX(degG);
  }

  //SetCoeff(G, 0, to_ZZ_p(181));
  //SetCoeff(G, 1, to_ZZ_p(451));
  //SetCoeff(G, 2, to_ZZ_p(528));
  //SetCoeff(G, 3, to_ZZ_p(611));
  // Make sure G is squarefree.
  while(!IsOne(GCD(G, diff(G)))){
    G = ZZ_pX(degG, 1) + random_ZZ_pX(degG);
  }

  // We set H(x) = x if deg(H) = 1 since that simplifies the arithmetic.
  if(degH == 0){
    H = ZZ_pX(0,1);
    D = G;
  }
  else if(degH == 1){
    SetX(H);
    D = LeftShift(G, 2);
    // Make sure x does not divide G.
    if(!IsOne(GCD(G,H))){
      while( !IsOne(GCD(G, diff(G))) && !IsOne(GCD(G,H)) ){
	G = ZZ_pX(degG, 1) + random_ZZ_pX(degG);
      }
    }
  }
  else{
    H = ZZ_pX(degH, 1) + random_ZZ_pX(degH);
    if(irred){
      do{
	while(!DetIrredTest(H))
	  H = ZZ_pX(degH, 1) + random_ZZ_pX(degH);
      }while(G == H);
    }
    // Make sure G and H are relatively prime and H is squarefree.
    while(!IsOne(GCD(G,H)) || !IsOne(GCD(H, diff(H)))  ){
      H = ZZ_pX(degH, 1) + random_ZZ_pX(degH);
    }
    D = G*sqr(H);
  }
  
  degD = deg(D);

  if(degD > 20){
    cout<<"Error: The degree of D(x) must be less than 21.\n\n";
    return;
  }

  // It may be that G and H are irreducible anyway.
  if(!irred){
    if(degH >=2){
      if(DetIrredTest(H) && DetIrredTest(G))
	irred = 1;
    }
    else{
      if(DetIrredTest(G))
	irred = 1;
    }
  }

  // Find the numbers and degrees of the factors of G and H.
  for(i=0; i<degD; i++)
    factordegs[i]=0;

  if(irred){
    factordegs[degG-1]++;
    if(degH)
      factordegs[degH-1]++;
  }
  else{
    if(degH==1)
      factordegs[0]++;
    else{
      factors = berlekamp(inv(LeadCoeff(H))*H);
      for(i=0; i<factors.length(); i++){
	factordegs[deg(factors[i].a)-1]++;
      }    
    }
    factors = berlekamp(inv(LeadCoeff(G))*G);
    
    for(i=0; i<factors.length(); i++){
      factordegs[deg(factors[i].a)-1]++;
    } 
  }
  return;
}

int main(int argc, char *argv[]){
  ZZ factor;
  int i;
  ZZ /*sumh = ZZ(),*/ sumE1diff = ZZ(), sumE2diff = ZZ();
  RR sumRatio = RR();
  vec_pair_ZZ_pX_long factors;

  if(argc != 2){
    cout<<"Usage: ./everything inputfile > outputfile\n\n";
    return 0;
  }
  else {
    if(!getinput(argv[1])) {
      return(0);
    }
  }

  SetSeed( to_ZZ( time(NULL) ) );

  setConstants();

  if(rank != 0){
    cout<<"Procedures not yet implemented to work with Rank "<<rank<<" curves.\n";
    return(0);
  }

  cout<<"\nComputing statistics for the distribution of:"<<endl;
  cout<<"h, |h-E_1| vs. |h-E_2|, and |h-E_2|/L^2.\n\n";
  cout<<"Taking "<<samples<<" samples of genus "<<genus<<" cubic curves over F_"<<q<<" with lambda = "<<lambda<<".\n\n";

  cout<<"h \t |h-E_1| \t |h-E_2| \t |h-E_1|/L_1^2 \t |h-E_2|/L_2^2 \t |h-E_2|/L_3^2  \t E1 \t E2 \t L_1^2 \t L_2^2 \t L_3^2 \t E_1-L_1^2 \t E_2-L_2^2 \t E_2-L_3^2 \t E_1+L_1^2 \t E_2+L_2^2 \t E_2+L_3^2 \t E_1-alpha_1*L_1^2 \t E_2-alpha_2*L_2^2 \t E_2-alpha_3*L_3^2 \t E_1+alpha_1*L_1^2 \t E_2+alpha_2*L_2^2 \t E_2+alpha_3*L_3^2 \t Alpha Int. 1 \t Alpha Int. 2\t Alpha Int. 3 \t Phase 1 \t BS time 1 \t GS time 1 \t BS time 2 \t GS time 2 \t BS time 3 \t GS time 3 \t B steps 1 \t G steps 1 \t B steps 2 \t G steps 2 \t B steps 3 \t G steps 3 \t Total Steps 1 \t Total Steps 2 \t Total Steps 3 \t D(x)\n";

  //BuildIrred(G, degG);
  //H = ZZ_pX(0,1);
  //D = G;

  for(i=0; i<samples; i++){
    // Step 0: Get random polynomials G and H
      
    getPolynomial();
    //BuildRandomIrred(G,D);
    //D = G;

    //G = ZZ_pX(0,3)+ZZ_pX(1, 1)+ZZ_pX(2, 15)+ZZ_pX(3, 3)+ZZ_pX(4, 14)+ZZ_pX(5, 1);

    //D = G*H*H;
    //cout<<"G(x) = "<<G<<endl;

    // Step 1: Compute an approximation E of the class number h
    //         and an integer L such that |h - E| < L^2.

    set1 = time(NULL);
    approxh();
    time1 = time(NULL);

    // Step 2: Use extra information about h in (E-L^2, E+L^2),
    //         such as its distribution in the interval
    //         or h mod r for small primes r.
    
    // Find elements of small order.
    //if(rank == 0){
    factor = smallOrder0();
    //cout<<"factor = "<<factor<<endl;
    //}

    // Step 3: Find h in the interval (E-L^2, E+L^2) via
    //         A) Baby Step, Giant Step or
    //         B) Pollard's Kangaroo
    //if(rank == 0){

    //cout<<"Beginning baby steps 0."<<endl;
    bsgs0(factor,0);
    //cout<<"Beginning baby steps 1."<<endl;
    bsgs0(factor,1);
    //cout<<"Beginning baby steps 2."<<endl;
    bsgs0(factor,2);
    //cout<<"Completed baby steps."<<endl;
     
    /*}
      else if(rank == 1)
      bsgs1();
      else
      bsgs2();
    */

    //sumh += h;
    if(IsZero(h)){
      i--;
    }
    else{
      
      //sumE1diff += abs(h-E1);
      //sumE2diff += abs(h-E);
      //sumRatio  += to_RR(abs(h-E))/to_RR(L22);
      
      // Print the data.
      cout<<h<<"\t"<<abs(h-Es[0])<<"\t"<<abs(h-Es[1])<<"\t"<<to_RR(abs(h-Es[0]))/to_RR(L2[0])<<"\t"<<to_RR(abs(h-Es[1]))/to_RR(L2[1])<<"\t"<<to_RR(abs(h-Es[1]))/to_RR(L2[2])<<"\t"<<Es[0]<<"\t"<<Es[1]<<"\t"<<L2[0]<<"\t"<<L2[1]<<"\t"<<L2[2]<<"\t"<<Es[0]-L2[0]<<"\t"<<Es[1]-L2[1]<<"\t"<<Es[1]-L2[2]<<"\t"<<Es[0]+L2[0]<<"\t"<<Es[1]+L2[1]<<"\t"<<Es[1]+L2[2]<<"\t"<<Es[0]-CeilToZZ(alpha1*to_RR(L2[0]))<<"\t"<<Es[1]-CeilToZZ(alpha2*to_RR(L2[1]))<<"\t"<<Es[1]-CeilToZZ(alpha3*to_RR(L2[2]))<<"\t"<<Es[0]+TruncToZZ(alpha1*to_RR(L2[0]))<<"\t"<<Es[1]+TruncToZZ(alpha2*to_RR(L2[1]))<<"\t"<<Es[1]+TruncToZZ(alpha3*to_RR(L2[2]))<<"\t"<<2*TruncToZZ(alpha1*to_RR(L2[0]))+1<<"\t"<<2*TruncToZZ(alpha2*to_RR(L2[1]))+1<<"\t"<<2*TruncToZZ(alpha3*to_RR(L2[2]))+1<<"\t"<<time1-set1<<"\t"<<timeb[0]-setb[0]<<"\t"<<timeg[0]-setg[0]<<"\t"<<timeb[1]-setb[1]<<"\t"<<timeg[1]-setg[1]<<"\t"<<timeb[2]-setb[2]<<"\t"<<timeg[2]-setg[2]<<"\t"<<babysteps[0]<<"\t"<<giantsteps[0]<<"\t"<<babysteps[1]<<"\t"<<giantsteps[1]<<"\t"<<babysteps[2]<<"\t"<<giantsteps[2]<<"\t"<<babysteps[0]+giantsteps[0]<<"\t"<<babysteps[1]+giantsteps[1]<<"\t"<<babysteps[2]+giantsteps[2]<<"\t\t\t"<<D<<endl;
    }
  }

  //cout<<"Total running time: "<<time1-set<<" seconds.\n";
  //cout<<"Average time per run = "<<(time1-set)/samples<<" seconds.\n\n";

  //cout<<"Average difference of |h-E_1| = "<<sumE1diff/samples<<".\n";
  //cout<<"Average difference of |h-E_2| = "<<sumE2diff/samples<<".\n";
  //cout<<"Average  ratio  |h-E_1|/L_2^2 = "<<sumRatio/samples<<".\n\n";

}
