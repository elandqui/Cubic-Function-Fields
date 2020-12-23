/************************************************************/
/* kangaroo_mpi.cc                                          */
/*                                                          */ 
/* This is the parallelized kangaroo algorithm written for  */
/* use with MPI on a cluster of processors.                 */
/* The function terminates when the master kangaroo writes  */
/* a file containing the solution.                          */
/*                                                          */ 
/* Written by Eric Landquist                                */
/************************************************************/


#include "kangaroo_mpi.h"

NTL_CLIENT

/********************************************************************/
/* Finds an approximation E of h and a bound L such that |h-E|<L^2. */
/********************************************************************/

void approxh(){
  int s1, s2;
  int n, v, divlim; // v|n
  int a, i, j, l;
  RR B;  // Set logE2 = A(K) + B.
  // ZZ B1; // Used to find the trivial bound E1.
  RR C, S, T; // Temp variables.
  RR Q = to_RR(q);
  RR A;
  RR psi3; //, psi2, psi1;
  ZZ SvSum; // Used for psi3.
  ofstream sfile;
  char *filename;
  fstream roofile;
  bool closed = true;
  ZZ split1, split2;
  time_t set1, time1;
  int i1time;

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

  // First compute logE2 = A(K) + B + C. 
  // Begin with A(K).
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
  Epr = exp(A);

  /***********************************************/
  /* Precompute the S_v(j) for 1 <= v <= lambda. */
  /***********************************************/

  // This will be the longest running step 
  // before the kangaroo jumping portion.
 
  sfile.open("getS.roo", ios::out);
  sfile<<q<<"\n"<<G<<"\n"<<H<<"\n"<<m<<"\n"<<lambda<<endl;
  sfile.close();
  
  int go_signal = 1;

  for(i=1; i<=m; i++)
    MPI_Send(&go_signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
  
  // The master will get the degree < lambda polynomials
  // The slaves will get the degree = lambda polynomials.
  set1 = time(NULL);
  for(v=1; v<=lambda-1; v++){
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
    Epr*=(pow(power(Q,2*v)/(power(Q,2*v) + power(Q,v) + 1.0),to_RR(splitting[1]))*pow(power(Q,2*v)/(power(Q,2*v) -2.0*power(Q,v) + 1.0),to_RR(splitting[2])));
  }
  time1 = time(NULL);
  
  phase1time = time1-set1;

  // Gather up the degree lambda results from the slaves.
  splitting[0] = splitting[1] = splitting[2] = 0;
  MPI_Status rc;
  int block;
  for(i=1; i<=m; i++){
    
    MPI_Recv(&block, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &rc);
    
    std::stringstream ss;
    filename = (char *)malloc(sizeof(char)*15);
    ss << "block" << block << ".roo";
    ss >> filename;

    while(closed){ 
      fstream roofile;
      roofile.open(filename, ios::in);
      if(roofile.is_open()) closed = false;
      else{
        roofile.close();
        sleep(1);
      }
    }

    roofile>>split1;
    roofile>>split2;
    roofile>>i1time;

    roofile.close();
    closed = true;
    
    splitting[1]+=split1;
    splitting[2]+=split2;
    phase1time+=i1time;
  }
  
  v = lambda;

  // Get the S-values needed to compute psi3.
  SvCache[2*(lambda-1)] = 2*splitting[2]-splitting[1];
  Sv1(lambda, 3);

  /*************************/
  /* Approximate h with E. */
  /*************************/  

  Epr*=(pow(power(Q,2*v)/(power(Q,2*v) + power(Q,v) + 1.0), to_RR(splitting[1]))*pow(power(Q,2*v)/(power(Q,2*v) -2.0*power(Q,v) + 1.0),to_RR(splitting[2])));

  E = RoundToZZ(Epr);

  /************************************************/
  /* Approximate the square root of the error, L. */
  /************************************************/

  // The first step is estimating the first term of this:
  // S_{lambda+1}(1)/(q^{lambda+1}) + ... + S_{2*lambda}(1)/(q^{2*lambda})
  if(q%3 == 1){
    //C = power(Q,lambda);
    T = sqrt(Q);

    i = lambda+1;
    n = i+1;

    /*************************************/
    /* These are the bounds: E3 and L3^2 */
    /*************************************/

    // E3 = E2.

    // This is a little bit of a speed-up. Instead of dividing odd i by 3,
    // and even i by 2, we divide i by its smallest prime factor.
    // There isn't a speed-up until i=5 is reached, i.e. unless lambda > 2.
    
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

    L32 = CeilToZZ(Epr*expm1(psi3)+1.0/2.0);

    // I'm going with this guy to speed things up.
    L = CeilToZZ(sqrt(alpha*to_RR(L32)));

    // Previously, there was a correction if the new interval
    // went out of the bounds of the Hasse interval.
    // That approach changed E, so we won't do that anymore.
    
    //cout<<"Hasse Interval: ["<<CeilToZZ(power(T - 1, 2*genus))<<", "<<FloorToZZ(power(T + 1, 2*genus))<<"].\n";
    //cout<<"New interval 3: ["<<E-L32<<", "<<E+L32<<"].\n\n";
    //cout<<"Estimate of the class number   E = "<<E<<".\n";
    //cout<<"Our bounds on the error:   L_3^2 = "<<L32<<".\n";
    //cout<<"Setting L = "<<L<<".\n\n";

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

/***********************************************************/
/* Get the splitting information for a set of polynomials. */
/***********************************************************/

void getSplitting(){

  fstream sfile;
  char *filename;
  bool closed = true;
  int block = 2*roonum + (rootype=='t' ? 0 : 1);
  ZZ split1, split2;
  vec_pair_ZZ_pX_long factors;
  int i;
  ZZ c, j, begin, end;
  ZZ_pX P;
  time_t set1, time1;
  
  begin = ((block*(q-1))/m)+1;
  if(block == 0)
    begin--;
  end = ((block+1)*(q-1))/m;

  for(i=0; i<21; i++) 
    factordegs[i]=0;
  
  // Get the factors of D(x).
  factors = berlekamp(inv(LeadCoeff(G))*G);
  for(i=0; i<factors.length(); i++){
    factordegs[deg(factors[i].a)-1]++;
  }    
  factors = berlekamp(inv(LeadCoeff(H))*H);
  for(i=0; i<factors.length(); i++){
    factordegs[deg(factors[i].a)-1]++;
  }    
  split1 = split2 = ZZ();
  P = ZZ_pX(lambda, 1);

  set1 = time(NULL);

  // Loop through every monic irreducible degree v polynomial.
  // In the case that deg(P) = v = 1, every polynomial is irred.
  if(lambda == 1){
    // If there are some ramified primes of this degree, 
    // We need to check those.
    if(factordegs[0]){
      for(j=begin; j<=end; j++){
	// Get the next polynomial.
	SetCoeff(P, 0, to_ZZ_p(j));
	if(IsZero(D%P)); // This test can be sped up.
	else if(chi(P) == 1){
	  // This is the complete splitting case.
	  split2++;
	}
	else{
	  // This is the inert case.
	  split1++;
	}
      }
    }
    else{
      for(j=begin; j<=end; j++){
	// Get the next polynomial.
	SetCoeff(P, 0, to_ZZ_p(j));
	if(chi(P) == 1){
	  // This is the complete splitting case.
	  split2++;
	}
	else{
	 // This is the inert case.
	  split1++;
	}
      }
    }
  }
  // We use a nice trick to run through every irreducible 
  // degree 2 polynomial.
  else if (lambda == 2){
    //P = ZZ_pX(v, 1);
    // If there are some ramified primes of this degree, 
    // We need to check those.
    if(factordegs[1]){
      for(j = begin; j <= end; j++){
	// If x^2 - j is irreducible, cycle through each (x+c)^2 - j.
	if(Jacobi(j, q) == -1){
	  SetCoeff(P, 0, to_ZZ_p(-j));
	  for(c = ZZ(); c < q; c++){
	    if(IsZero(D%P)); // This test can be sped up.
	    else if(chi(P) == 1){
	      // This is the complete splitting case.
	      split2++;
	    }
	    else{
	      // This is the inert case.
	      split1++;
	    }
	    // Get the next polynomial.
	    SetCoeff(P, 0, to_ZZ_p(rep(coeff(P,0))+(2*c+1)));
	    SetCoeff(P, 1, to_ZZ_p(rep(coeff(P,1))+2));
	  }
	}
      }
    }
    else{
      for(j = begin; j <= end; j++){
	// If x^2 - j is irreducible, cycle through each (x+c)^2 - j.
	if(Jacobi(j, q) == -1){
	  SetCoeff(P, 0, to_ZZ_p(-j));
	  for(c = ZZ(); c < q; c++){
	    if(chi(P) == 1){
	      // This is the complete splitting case.
	      split2++;
	    }
	    else{
	      // This is the inert case.
	      split1++;
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
    SetCoeff(P, lambda-1, to_ZZ_p(begin));
    // Here we have deg(P) > 2, so we test for irreducibility.
    if(factordegs[lambda-1]){
      while(IsOne(LeadCoeff(P)) && (rep(coeff(P, lambda-1)) < end)){
	// Formerly if(ProbIrredTest(P, iter=1)){
	if(DetIrredTest(P)){
	  if(IsZero(D%P));
	  else if(chi(P) == 1){
	    split2++;
	  }
	  else{
	    split1++;
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
      while(IsOne(LeadCoeff(P)) && (rep(coeff(P, lambda-1)) < end) ){
	
	if(DetIrredTest(P)){
	  if(chi(P) == 1){
	    split2++;
	  }
	  else{
	    split1++;
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
  
  time1 = time(NULL);

  /**********************************/
  /* Write the results to the file. */
  /**********************************/
  int i1time = time1-set1;

  std::stringstream ss;
  filename = (char *)malloc(sizeof(char)*15);
  ss << "block" << block << ".roo";
  ss >> filename;

  sfile.open(filename, ios::out);
  sfile<<split1<<"\n"<<split2<<"\n"<<i1time;
  sfile.close();
   
  /*******************************************************/
  /* Notify the master kangaroo that this block is done. */
  /*******************************************************/

  MPI_Send(&block, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

/***********************************************/
/* Initialize the kangaroo jumps, unit rank 0. */
/***********************************************/

void getRooJumps0(){
  int i, found = 0;
 
  ZZ avgjump;        // The average of the jumps.
  ZZ randjump;       // The upper bound on the random jump distance.
  ZZ start;          // The starting position for t_0.
  ZZ offset = ZZ();  // Separate the kangaroos by offset units.
  cubic_ideal g;     // The generator
  cubic_ideal temp;
  cubic_ideal jumps[64]; 
  ZZ jumpdistance[64];
  ZZ jumpsum;
  
  /******************************************************/
  /* Set the average jump distance, the jump distances, */
  /* and the starting position for tame kangaroos.      */
  /******************************************************/
  RoundToZZ(avgjump, to_RR(m*L)*sqrt(to_RR(l))/2);
  avgjump-=(avgjump%l);
  randjump = 2*avgjump/l;
  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = l*(RandomBnd(randjump)+1);
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*avgjump) || (jumpsum > (i+1)*avgjump) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = l*(RandomBnd(randjump)+1);
      jumpsum += jumpdistance[i];
    }
    jumpdistance[63] -= jumpdistance[i];
  }
  
  jumpdistance[63] += 64*avgjump;
  
  start = E - E%l + a;

  if(m > 2){
    offset = avgjump/((m/2)-1);
    offset -= (offset % l);
  }

  // Make sure that the offset isn't a jump distance.
  i=0;
  while(i<64){
    if(offset == jumpdistance[i]){
      offset += l;
      i = 0;
    }
    else
      i++;
  }

  // Set the jumps. First generate a random generator.
  random(g);
  for(i=0; i<64; i++)
    jumps[i] = g^jumpdistance[i];

  //cout<<avgjump<<endl;

  // Set the number of bits for a distinguished ideal.
  distbits = (NumBits(avgjump)-1)/2 + 2;// (NumBits(avgjump)-1)/2;
  distbits += (distbits%2);

  /****************************************/
  /* Print start and distbits to ini.roo. */
  /****************************************/
  ofstream roofile;

  roofile.open("ini.roo", ios::out);
  //roofile<<q<<"\n"<<G<<"\n"<<H<<"\n"<<m<<endl;
  roofile<<start<<"\n"<<offset<<"\n"<<distbits<<endl;

  /************************************/
  /* Print jump distances to roo.ini. */
  /************************************/

  for(i=0; i<64; i++)
    roofile<<jumpdistance[i]<<endl;

  roofile<<g.s<<"\n"<<g.s1<<"\n"<<g.s2<<"\n"<<g.u<<"\n"<<g.v<<"\n"<<g.w<<endl;

  for(i=0; i<64; i++)
    roofile<<jumps[i].s<<"\n"<<jumps[i].s1<<"\n"<<jumps[i].s2<<"\n"<<jumps[i].u<<"\n"<<jumps[i].v<<"\n"<<jumps[i].w<<endl;

  roofile.close();

  // Write the go file.
  //ofstream gofile;
  //gofile.open("go.roo",ios::out);
  //gofile<<"Boing! Boing! Boing!"<<endl;
  //gofile.close();

}

/*********************************************/
/* Get the kangaroo jumps: unit rank 1 case. */
/*********************************************/

void getRooJumps1(){
  int i, found = 0;
 
  ZZ avgjump;        // The average of the jumps.
  ZZ randjump;       // The upper bound on the random jump distance.
  ZZ offset = ZZ();  // The spacing between kangaroos.
  infrastructure_ideal temp;
  cubic_ideal jumps[64]; 
  ZZ jumpdistance[64];
  ZZ jumpsum = ZZ();
  long headwind;
  // Get the headwind.
  headwind1(headwind);

  avgjump = CeilToZZ((double)m*sqrt(alpha*to_RR(U)*(2.0*tauratio-1.0)) - 2.0*(tauratio-1.0));
   
  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 2 + headwind);
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*avgjump) || (jumpsum > (i+1)*avgjump) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 2 + headwind);
      jumpsum += jumpdistance[i];
    }
    jumpdistance[63] -= jumpdistance[i];
  }
  
  jumpdistance[63] += (64*avgjump + headwind + 32);
  
  if(m > 2)
    offset = avgjump/((m/2)-1);

  // Set the jumps. 
  for(i=0; i<64; i++){
    below(jumpdistance[i], temp);
    jumpdistance[i] = temp.d0;
    inf_to_cubic(temp, jumps[i]);
  }

  // Set the number of bits for a distinguished ideal.
  distbits = (NumBits(avgjump)-1)/2;// (NumBits(avgjump)-1)/2;
  distbits += (distbits%2);

  /****************************************/
  /* Print start and distbits to roo.ini. */
  /****************************************/
  ofstream roofile;

  roofile.open("ini.roo", ios::out);
  //roofile<<q<<"\n"<<G<<"\n"<<H<<"\n"<<m<<endl;
  roofile<<E<<"\n"<<offset<<"\n"<<distbits<<endl;

  /************************************/
  /* Print jump distances to roo.ini. */
  /************************************/

  for(i=0; i<64; i++)
    roofile<<jumpdistance[i]<<endl;

  for(i=0; i<64; i++)
    roofile<<jumps[i].s<<"\n"<<jumps[i].s1<<"\n"<<jumps[i].s2<<"\n"<<jumps[i].u<<"\n"<<jumps[i].v<<"\n"<<jumps[i].w<<endl;

  roofile.close();
}

/*********************************************/
/* Get the kangaroo jumps: unit rank 2 case. */
/*********************************************/

void getRooJumps2(){
  int i, found = 0;
 
  ZZ avgjump;        // The average of the jumps.
  ZZ randjump;       // The upper bound on the random jump distance.
  ZZ offset = ZZ();  // The spacing between kangaroos.
  infrastructure_ideal temp;
  cubic_ideal jumps[64]; 
  ZZ jumpdistance[64];
  ZZ jumpdistance1[64];
  ZZ jumpsum = ZZ();
  long dist0, dist1, dist2;
  ZZ wind0, wind2; 

  if(tauratio <= to_double(l)){
    tauratio = 1.0;
  }
  else{
    tauratio /= to_double(l);
  }

  // Get the headwind.
  headwind2(dist0, dist1, dist2);
  wind0 = to_ZZ(dist0);
  wind2 = to_ZZ(dist2);

  avgjump = CeilToZZ((double)m*sqrt(alpha*to_RR(U*l)*(2.0*tauratio-1.0))/2.0 - to_RR(l)*(tauratio-1.0));
  randjump = (2*avgjump - (genus+1))/l;
  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = l*(RandomBnd(randjump) + genus + 1) + wind0;
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*avgjump) || (jumpsum > (i+1)*avgjump) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = l*(RandomBnd(randjump) + genus + 1) + wind0;
      jumpsum += jumpdistance[i];
    }
    jumpdistance[63] -= jumpdistance[i];
  }
  jumpdistance[63] += 64*(avgjump);
  jumpdistance[63] -= jumpdistance[63]%l;
  jumpdistance[63] += wind0;
  
  if(m > 2){
    offset = avgjump/((m/2)-1);
    offset -= (offset%l);
  }
  
  // Set the jumps. 
  for(i=0; i<64; i++){
    below(jumpdistance[i], wind2, temp);
    while((temp.d0%l != wind0) && (temp.d2 != wind2) ){
      jumpdistance[i]+=l;
      below(jumpdistance[i], wind2, temp);
    } 
    jumpdistance[i] = temp.d0;
    jumpdistance1[i] = temp.d1;
    inf_to_cubic(temp, jumps[i]);
  }

  // Set the number of bits for a distinguished ideal.
  distbits = (NumBits(avgjump)-1)/2;// (NumBits(avgjump)-1)/2;
  distbits += (distbits%2);

  /****************************************/
  /* Print start and distbits to roo.ini. */
  /****************************************/
  ofstream roofile;
  E -= E%l;
  E += a;
  roofile.open("ini.roo", ios::out);
  //roofile<<q<<"\n"<<G<<"\n"<<H<<"\n"<<m<<endl;
  roofile<<E<<"\n"<<offset<<"\n"<<distbits<<endl;

  /************************************/
  /* Print jump distances to roo.ini. */
  /************************************/

  for(i=0; i<64; i++)
    roofile<<jumpdistance[i]<<"\n"<<jumpdistance1[i]<<endl;

  for(i=0; i<64; i++)
    roofile<<jumps[i].s<<"\n"<<jumps[i].s1<<"\n"<<jumps[i].s2<<"\n"<<jumps[i].u<<"\n"<<jumps[i].v<<"\n"<<jumps[i].w<<endl;

  roofile.close();
}

/***************************************/
/* Reads the distinguished point file. */
/***************************************/

int readDPfile(){
  ifstream dpfile, lastdpfile;
  long numdps = 0;
  //bool closed = true;
  int i=0;
  int matchloc;
  //ZZX input1;
  //ZZ input2;
  char *filename, *lastdpfilename;
  ZZ d1, d2; // dummy variables
  
  type = (char*)malloc(sizeof(char)*size);
  nums = (int*)malloc(sizeof(int)*size);

  // Let's get the contents.
  for(int k=0; k<mt+mw; k++){
    std::stringstream ss, ss2;
    filename = (char *)malloc(sizeof(char)*17);
    lastdpfilename = (char *)malloc(sizeof(char)*18);

    if(k<mt){
      ss << "all_dp_t" << k << ".roo";
      ss >> filename;
      ss2 << "last_dp_t"<< k << ".roo";
      ss2 >> lastdpfilename;
    }
    else{
      ss << "all_dp_w" << k-mt << ".roo";
      ss >> filename;
      ss2 << "last_dp_w"<< k-mt << ".roo";
      ss2 >> lastdpfilename;
    }
    //closed = true;
    //while(closed){
    dpfile.open(filename, ios::in);
    if( dpfile.is_open() ){ 
      
      // Check to see how many distinguished points
      // this kangaroo has found.
      
      lastdpfile.open(lastdpfilename, ios::in);
      if( lastdpfile.is_open() )
	lastdpfile >> numdps;
      else
	numdps = 0;
      lastdpfile.close();
      
      // Let's get the contents.
      
      //while(!dpfile.eof()){
      for(j=0; j<numdps; j++){ 
	dpfile >> dp[i];
	dpfile >> dists[i];
	if(rank == 2){
	  dpfile >> d1;
	  dpfile >> d2;
	}
	if(k<mt){ 
	  type[i] = 't';
	  roonum[i] = k;
	}
	else{ 
	  type[i] = 'w';
	  roonum[i] = k-mt;
	}
	i++;
	
	// Need bigger arrays!
	if(size == i){
	  size+=100;
	  dp.SetLength(size);
	  dists.SetLength(size);
	  type = (char*)realloc(type, sizeof(char)*size);
	  roonum = (int*)realloc(roonum, sizeof(int)*size);
	}
      }  
      //closed = false;
    }
    dpfile.close();
    //else  {
    //	dpfile.close();
    //	sleep(1);
    //}  
    //}
  }
    
  // Sort the distinguished point arrays.
  sort(0, i-1);

  // Search for matches.
  return search(i, matchloc);
}

/***************************************/
/* Sort the distinguished point array. */
/***************************************/

void sort(int lo, int hi){

  //  lo is the lower index, hi is the upper index
  //  of the region of array a that is to be sorted
  int i=lo, j=hi;
  ZZX swap1;
  ZZ swap2;
  char swap3;
  int swap4;
  ZZ x = coeff(dp[(lo+hi)/2], 0);

  //  partition
  do  {    
    while ( coeff(dp[i],0) < x ) i++; 
    while ( coeff(dp[j],0) > x ) j--;
    if (i<=j){
      swap1 = dp[i]; dp[i] = dp[j]; dp[j] = swap1;
      swap2 = dists[i]; dists[i] = dists[j]; dists[j] = swap2;
      swap3 = type[i]; type[i] = type[j]; type[j] = swap3;
      swap4 = nums[i]; nums[i] = nums[j]; nums[j] = swap4;
      i++; j--;
    }
  } while (i<=j);

  //  recursion
  if (lo<j) sort(lo, j);
  if (i<hi) sort(i, hi);
}


/***********************/
/* Search for matches. */
/***********************/

int search(int len, int &matchloc){
  int i;

  for(i=1; i<len; i++){
    if(dp[i-1]==dp[i]){
      // Two kangaroos of the same type collided. Reassign one
      if(type[i-1]==type[i]){
        reassign(type[i-1], nums[i-1]);
      }
      // This is a match!
      else{
	if(rank == 0){
	  order = abs(dists[i-1]-dists[i]);
	  if(check(order)){
	    matchloc = i;
	    return(1);
	  }
	}
	else if (rank == 1){
	  order = abs(dists[i-1]-dists[i])/2; 
	  if(check_inf(order)){
	    matchloc = i;
	    return(1);
	  }
	}
	else { // rank == 2.
	  order = abs(dists[i-1]-dists[i]); 
	  if(check_inf2(order)){
	    matchloc = i;
	    return(1);
	  }
	}
      }
    }
  }

  // If it gets here, there is no match.
  return(0);
}

/*************************/
/* Reassigns a kangaroo. */
/*************************/

void reassign(char tw, int num){
  char *rooname;
  fstream roofile;

  std::stringstream ss;
  rooname = (char *)malloc(18*sizeof(char));

  ss << "reassign_" << tw << num << ".roo";
  ss >> rooname;
  // First check to see if this kangaroo has been reassigned.
  roofile.open(rooname,ios::in);
  if( roofile.is_open() )  {
    // If the reassign file exists, then this roo has been reassigned.
    roofile.close();
    return;
  }
  roofile.close();
  roofile.open(rooname,ios::out);
  // Write the name of the new kangaroo.
  if(tw=='t'){
    roofile<<mt<<endl;
    mt++;
  }
  else{
    roofile<<mw<<endl;
    mw++;
  }
  
  // Write the shift from the current position.
  if(rank == 0)
    roofile<<l<<endl;
  else if (rank == 1)
    roofile<<3*tauratio<<endl;
  else{
    long wind0, wind1, wind2;
    headwind2(wind0, wind1, wind2);
    roofile<<RoundToZZ(3.0*tauratio*to_RR(l))+wind0<<endl;
  }
  roofile.close();
}

/*******************/
/* Make the jumps. */
/*******************/

void kangaroo0(){
  int i, found = 0;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  cubic_ideal roo;     // The kangaroo, tame or wild.
  cubic_ideal temp;
  cubic_ideal g;
  cubic_ideal jumps[64];
  int distbits2 = distbits/2;
  ZZ distance = ZZ();
  ZZ jumptotal = ZZ(); // The number of jumps for this kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  fstream dpfile;

  /***************************************/
  /* Initialize the generator and jumps. */
  /***************************************/

  g.s = to_ZZ_pX(gz[0]);
  g.s1 = to_ZZ_pX(gz[1]);
  g.s2 = to_ZZ_pX(gz[2]);
  g.u = to_ZZ_pX(gz[3]);
  g.v = to_ZZ_pX(gz[4]);
  g.w = to_ZZ_pX(gz[5]);

  for(i=0; i<64; i++){
    jumps[i].s = to_ZZ_pX(jumpsz[6*i]);
    jumps[i].s1 = to_ZZ_pX(jumpsz[6*i+1]);
    jumps[i].s2 = to_ZZ_pX(jumpsz[6*i+2]);
    jumps[i].u = to_ZZ_pX(jumpsz[6*i+3]);
    jumps[i].v = to_ZZ_pX(jumpsz[6*i+4]);
    jumps[i].w = to_ZZ_pX(jumpsz[6*i+5]);
  }
  
  /****************************/
  /* Initialize the kangaroo. */
  /****************************/
  
  // If we're picking up where we left off, then get the initialization data.
  dpfile.open(lastdpname, ios::in);
  if( dpfile.is_open() )  {
    // Read information in.
    dpfile >> distpoints;
    dpfile >> jumptotal;
    dpfile >> roo.s;
    dpfile >> roo.s1;
    dpfile >> roo.s2;
    dpfile >> roo.u;
    dpfile >> roo.v;
    dpfile >> roo.w;
    dpfile >> distance;
    dpfile.close();
  }
  else{
    dpfile.close();
    // This is a tame kangaroo.
    if(rootype=='t'){
      // Set the tame kangaroo.
      distance = start + offset*roonum;
    }
    // This is a wild kangaroo.
    else{
      distance = offset*roonum;
    }

    roo = g^distance;
  }

  /************************************************/
  /* This is the heart of the kangaroo algorithm. */
  /************************************************/

  while(!found){

    /***********************************/
    /* Check for distinguished points. */
    /***********************************/

    roocoeff = rep(coeff(roo.s,0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(roo.s,1));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
        distpoints++;
        found = distpoint(roo, distance, distpoints, jumptotal); 
	//cout<<"New DP: "<<distpoints<<" "<<found<<endl;
      }
    }

    /******************/
    /* Make the jump. */
    /******************/

    i = vmap(roo);
    temp = cubic_ideal(roo);
    roo = temp*jumps[i];
    distance += jumpdistance[i];
    jumptotal++;

  }

}

/*************************************/
/* Make the jumps. Unit rank 1 case. */
/*************************************/

void kangaroo1(){
  int i, found = 0;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  infrastructure_ideal roo;     // The kangaroo, tame or wild.
  infrastructure_ideal temp;
  cubic_ideal jumps[64];
  int distbits2 = distbits/2;
  ZZ distance = ZZ();
  //ZZ jumptotal = ZZ(); // The number of jumps for this kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  ZZ bsjumptotal = ZZ(), gsjumptotal = ZZ();
  fstream dpfile;
  ZZ Stau = RoundToZZ(to_RR(q)/tau3());

  /*************************/
  /* Initialize the jumps. */
  /*************************/

  for(i=0; i<64; i++){
    jumps[i].s  = to_ZZ_pX(jumpsz[6*i]);
    jumps[i].s1 = to_ZZ_pX(jumpsz[6*i+1]);
    jumps[i].s2 = to_ZZ_pX(jumpsz[6*i+2]);
    jumps[i].u  = to_ZZ_pX(jumpsz[6*i+3]);
    jumps[i].v  = to_ZZ_pX(jumpsz[6*i+4]);
    jumps[i].w  = to_ZZ_pX(jumpsz[6*i+5]);
  }
  
  /****************************/
  /* Initialize the kangaroo. */
  /****************************/
  
  // If we're picking up where we left off, then get the initialization data.
  dpfile.open(lastdpname, ios::in);
  if( dpfile.is_open() )  {
    // Read information in.
    dpfile >> distpoints;
    dpfile >> gsjumptotal;
    dpfile >> bsjumptotal;
    dpfile >> roo.mu0;
    dpfile >> roo.mu1;
    dpfile >> roo.mu2;
    dpfile >> roo.nu0;
    dpfile >> roo.nu1;
    dpfile >> roo.nu2;
    dpfile >> roo.d;
    dpfile >> roo.d0;
    dpfile.close();
  }
  else{
    dpfile.close();
    // This is a tame kangaroo.
    if(rootype=='t'){
      // Set the tame kangaroo.
      distance = 2*(start + offset*roonum);
    }
    // This is a wild kangaroo.
    else{
      distance = 2*offset*roonum;
    }
    below(distance, roo);
  }
  
  while(!found){
    /************************************************/
    /* Make baby steps until L(roo)(0) < [q/tau_3]. */
    /************************************************/
    while(rep(coeff(roo.d, 0)) >= Stau){
      baby_step_r1(roo, temp);
      roo = infrastructure_ideal(temp);

      bsjumptotal++;
    }
    
    /***********************************/
    /* Check for distinguished points. */
    /***********************************/
    // I'll only check for DPs in the set S_tau.
    roocoeff = rep(coeff(roo.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(roo.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	distpoints++;
	found = distpoint(roo, distpoints, bsjumptotal, gsjumptotal);//, false);
      }
    }
    
    /******************/
    /* Make the jump. */
    /******************/

    i = vmap(roo);
    temp = infrastructure_ideal(roo);
    giant_step_r1(roo, temp, jumps[i], jumpdistance[i]);
    //roo = temp*jumps[i];
    gsjumptotal++;
  }
}

/*************************************/
/* Make the jumps. Unit rank 2 case. */
/*************************************/

void kangaroo2(){
  int i, j, found = 0, count = 0;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  infrastructure_ideal roo;     // The kangaroo, tame or wild.
  infrastructure_ideal temp;
  cubic_ideal jumps[64];
  int distbits2 = distbits/2;
  ZZ distance = Zero;
  //ZZ jumptotal = ZZ(); // The number of jumps for this kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  ZZ bsjumptotal = Zero, gsjumptotal = Zero;
  fstream dpfile;
  ZZ Stau;
  long dist0, dist1, dist2;
  ZZ wind0, wind2; 
  ZZ tempd1;
  ZZ eclass;
  //bool checkLattice = true;
  //ZZ icn = One;
  // Make sure the roo hops in the correct equivalence class.
  if(rootype == 't')
    eclass = a;
  else
    eclass = Zero;

  if(tauratio <= to_double(l)){
    tauratio = 1.0;
  }
  else{
    tauratio /= to_double(l);
  }
  Stau = RoundToZZ(to_RR(q)/tauratio);

  headwind2(dist0, dist1, dist2);
  wind0 = to_ZZ(dist0);
  wind2 = to_ZZ(dist2);

  /*************************/
  /* Initialize the jumps. */
  /*************************/

  for(i=0; i<64; i++){
    jumps[i].s  = to_ZZ_pX(jumpsz[6*i]);
    jumps[i].s1 = to_ZZ_pX(jumpsz[6*i+1]);
    jumps[i].s2 = to_ZZ_pX(jumpsz[6*i+2]);
    jumps[i].u  = to_ZZ_pX(jumpsz[6*i+3]);
    jumps[i].v  = to_ZZ_pX(jumpsz[6*i+4]);
    jumps[i].w  = to_ZZ_pX(jumpsz[6*i+5]);
  }
  
  /****************************/
  /* Initialize the kangaroo. */
  /****************************/
  
  // If we're picking up where we left off, then get the initialization data.
  dpfile.open(lastdpname, ios::in);
  if( dpfile.is_open() )  {
    // Read information in.
    dpfile >> distpoints;
    dpfile >> gsjumptotal;
    dpfile >> bsjumptotal;
    dpfile >> roo.mu0;
    dpfile >> roo.mu1;
    dpfile >> roo.mu2;
    dpfile >> roo.nu0;
    dpfile >> roo.nu1;
    dpfile >> roo.nu2;
    dpfile >> roo.d;
    dpfile >> roo.d0;
    dpfile >> roo.d1;
    dpfile >> roo.d2;
    dpfile.close();
  }
  else{
    dpfile.close();
    // This is a tame kangaroo.
    if(rootype=='t'){
      // Set the tame kangaroo.
      distance = start + offset*roonum;
    }
    // This is a wild kangaroo.
    else{
      distance = offset*roonum;
    }
    below(distance, Zero, roo);
  }
  
  while(!found){
    /************************************************/
    /* Make baby steps until L(roo)(0) < [q/tau_3]. */
    /************************************************/
    while( (rep(coeff(roo.d, 0)) >= Stau) && (roo.d0%l == eclass) ){
      baby_step_0_r2(roo, temp);
      roo = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to 0 as possible.
      if(roo.d2 < 0){
	reduce_basis2(roo);
	while(roo.d2 <= 0){
	  temp = roo;
	  baby_step_2_r2(temp, roo);
	  count++;
	}
	reduce_basis(temp);
	for(j=0; j<count-1; j++){
	  roo = temp;
	  baby_step_0_r2(roo, temp);
	}
	roo = temp;
	count=0;
      }

      bsjumptotal++;
    }
    
    /***********************************/
    /* Check for distinguished points. */
    /***********************************/
    // I'll only check for DPs in the set S_tau.
    roocoeff = rep(coeff(roo.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(roo.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	distpoints++;	
	found = distpoint(roo, distpoints, bsjumptotal, gsjumptotal);
      }
    }
    
    /******************/
    /* Make the jump. */
    /******************/

    i = vmap(roo);
    temp = infrastructure_ideal(roo);
    giant_step_r2(roo, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    
    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(roo.d2 < 0){
      reduce_basis2(roo);
      while(roo.d2 <= 0){
	temp=roo;
	baby_step_2_r2(temp, roo);
	count++;
      }
      reduce_basis(temp);
      for(j=0; j<count-1; j++){
	roo = temp;
	baby_step_0_r2(roo, temp);
      }
      roo = temp;
      count=0;
    }
    else if(roo.d2 > 0){
      reduce_basis1(roo);
      while(roo.d2 >= 0){
	temp=roo;
	baby_step_1_r2(temp, roo);
      }
      reduce_basis(temp);
      roo = temp;
    }
    
    // Get the jumps back in the correct residue classes.

    while(roo.d0%l != eclass){
      temp = roo;
      baby_step_0_r2(temp, roo);
    }

    gsjumptotal++;
  }
}

/****************************************************/
/* Functions mapping a kangaroo to {0, 1, ..., 63}. */
/****************************************************/

inline int vmap(cubic_ideal &A){
  return(trunc_long(rep(coeff(A.u,0)), 6));
  //return(rep(coeff(A.u,0))%50);
}

inline int vmap(infrastructure_ideal &A){
  return(trunc_long(rep(coeff(A.mu1, 0)), 6));
}

/**********************************/
/* Process a distinguished point. */
/**********************************/

// Write the distinguished point to a file.
// Return 0 if a solution has not been found, and 1 if a solution has been found.
int distpoint(cubic_ideal &roo, ZZ &distance, long dps, ZZ &jumps){

  ofstream dpfile1, dpfile2;
  fstream roofile;
  fstream solnfile;
  bool closed = true;

  // Check to see if a solution has been found.

  solnfile.open("solution.txt",ios::in);
  if( solnfile.is_open() )  {
    // Solution file exists. Return 1 and end the jumps.
    solnfile.close();
    return(1);
  }
  solnfile.close();

  // If this kangaroo is being reassigned because it collided with another
  // kangaroo of the same type, then we read in the new information.

  std::stringstream ss;
  char *filename = (char *)malloc(18*sizeof(char));

  ss << "reassign_" << rootype << roonum << ".roo";
  ss >> filename;
  // First check to see if this kangaroo has been reassigned.
  roofile.open(filename,ios::in);
  if( roofile.is_open() )  {
    ZZ shift;
    roofile>>roonum;
    roofile>>shift;
    // If the reassign file exists, then this roo has been reassigned.
    roofile.close();
    distance+=shift;
    cubic_ideal g, shiftjump, temp;
    g.s = to_ZZ_pX(gz[0]);
    g.s1 = to_ZZ_pX(gz[1]);
    g.s2 = to_ZZ_pX(gz[2]);
    g.u = to_ZZ_pX(gz[3]);
    g.v = to_ZZ_pX(gz[4]);
    g.w = to_ZZ_pX(gz[5]);
    shiftjump = g^shift;
    temp = cubic_ideal(roo);
    roo = temp*shiftjump;
    
    jumps++;

    // Create the distinguished point filenames for this kangaroo.
    std::stringstream ss;
    lastdpname = (char *)malloc(17*sizeof(char));
    
    ss << "last_dp_" << rootype << roonum << ".roo";
    ss >> lastdpname;
    
    std::stringstream ss2;
    alldpname = (char *)malloc(16*sizeof(char));
    
    ss2 << "all_dp_" << rootype << roonum << ".roo";
    ss2 >> alldpname;
    return(0);
  }
  roofile.close();

  // Open the distinguished point file.
    
  while(closed){
    dpfile1.open(alldpname, ios::app);
    if( dpfile1.is_open() ) closed = false;
    else  {
      dpfile1.close();
      sleep(1);
    }  
  }

  // Write the distinguished points and distance to the alldp file.
  // Structure of dpoints.roo:
  // roo.s [return] distance [return] rootype [return] roonum [return]
  dpfile1<<"\n"<<roo.s<<"\n"<<distance;
 
  dpfile1.close();
 
  // Now write the distinguished point data to a file to read in case
  // the processor shuts off and we need to pick up where we left off.

  dpfile2.open(lastdpname, ios::out);
  dpfile2<<dps<<"\n"<<jumps<<"\n"<<roo.s<<"\n"<<roo.s1<<"\n"<<roo.s2<<"\n"<<roo.u<<"\n"<<roo.v<<"\n"<<roo.w<<"\n"<<distance<<endl;
  dpfile2.close();

  return(0);
}

/**********************************************************/
/* Process a distinguished point. Unit rank 1 and 2 case. */
/**********************************************************/

// Write the distinguished point to a file.
// Return 0 if a solution has not been found, and 1 if a solution has been found.
int distpoint(infrastructure_ideal &roo, long dps, ZZ &bsjumps, ZZ &gsjumps){//, bool &checkLattice){
  int j;
  ofstream dpfile1, dpfile2;
  fstream roofile;
  fstream solnfile;
  //fstream latticefile;
  bool closed = true;

  // Check to see if a solution has been found.
  solnfile.open("solution.txt",ios::in);
  if( solnfile.is_open() )  {
    // Solution file exists. Return 1 and end the jumps.
    solnfile.close();
    return(1);
  }
  solnfile.close();

  // In the unit rank 2 case, check to see if we need to
  // begin scanning the lattice for the regulator.
  /*
  if(checkLattice){
    latticefile.open("latticeini.roo",ios::in);
    if( latticefile.is_open() )  {
      // Lattice file exists. Dig out the data and begin jumping
      // with the new kangaroos.
      latticefile.close();
      checkLattice = false;
      return(0);
    }
    latticefile.close();
  }
  */
  // If this kangaroo is being reassigned because it collided with another
  // kangaroo of the same type, then we read in the new information.
  std::stringstream ss;
  char *filename = (char *)malloc(18*sizeof(char));

  ss << "reassign_" << rootype << roonum << ".roo";
  ss >> filename;

  // First check to see if this kangaroo has been reassigned.
  roofile.open(filename,ios::in);
  if( roofile.is_open() )  {
    ZZ shift;
    roofile>>roonum;
    roofile>>shift;
    // If the reassign file exists, then this roo has been reassigned.
    roofile.close();
    infrastructure_ideal shiftjump, temp;
    if(rank == 1)
      below(shift, shiftjump);
    else{
      long wind0, wind1, wind2;
      headwind2(wind0, wind1, wind2);
      ZZ wind = to_ZZ(wind2);
      below(shift, wind, shiftjump);
    }
    temp = infrastructure_ideal(roo);
    roo = temp*shiftjump;
    if(rank == 2){
      int count = 0;
      if(roo.d2 < 0){
	reduce_basis2(roo);
	while(roo.d2 <= 0){
	  temp = roo;
	  baby_step_2_r2(temp, roo);
	  count++;
	}
	reduce_basis(temp);
	for(j=0; j<count-1; j++){
	  roo = temp;
	  baby_step_0_r2(roo, temp);
	}
	roo = temp;
	count=0;
      }
    }
    while(roo.d0%l != eclass){
      temp = roo;
      baby_step_0_r2(temp, roo);
    }
    gsjumps++;

    // Create the distinguished point filenames for this kangaroo.
    std::stringstream ss;
    lastdpname = (char *)malloc(17*sizeof(char));
    
    ss << "last_dp_" << rootype << roonum << ".roo";
    ss >> lastdpname;
    
    std::stringstream ss2;
    alldpname = (char *)malloc(16*sizeof(char));
    
    ss2 << "all_dp_" << rootype << roonum << ".roo";
    ss2 >> alldpname;
    
    return(0);
  }
  roofile.close();

  // Open the distinguished point file.
  dpfile1.open(alldpname, ios::app);
  if( dpfile1.is_open() ){
    
    // Write the distinguished points and distance to the dpoints.roo file.
    // Structure of dpoints.roo:
    // roo.s [return] distance [return] rootype [return] roonum [return]
    if(rank == 1)
      dpfile1<<"\n"<<roo.d<<"\n"<<roo.d0;
    else
      dpfile1<<"\n"<<roo.d<<"\n"<<roo.d0<<"\n"<<roo.d1<<"\n"<<roo.d2;
    
    dpfile1.close();
    closed = false;
  }
  else  {
    dpfile1.close();
    sleep(1);
  }  
 
  // Now write the distinguished point data to its own file to read in case
  // the processor shuts off and we need to pick up where we left off.

  dpfile2.open(lastdpname, ios::out);
  if(rank == 1)
    dpfile2<<dps<<"\n"<<gsjumps<<"\n"<<bsjumps<<"\n"<<roo.mu0<<"\n"<<roo.mu1<<"\n"<<roo.mu2<<"\n"<<roo.nu0<<"\n"<<roo.nu1<<"\n"<<roo.nu2<<"\n"<<roo.d<<"\n"<<roo.d0<<endl;
  else
    dpfile2<<dps<<"\n"<<gsjumps<<"\n"<<bsjumps<<"\n"<<roo.mu0<<"\n"<<roo.mu1<<"\n"<<roo.mu2<<"\n"<<roo.nu0<<"\n"<<roo.nu1<<"\n"<<roo.nu2<<"\n"<<roo.d<<"\n"<<roo.d0<<"\n"<<roo.d1<<"\n"<<roo.d2<<endl;
  dpfile2.close();

  return(0);
}

/**************************************************/
/* Reads initialization information from ini.roo. */
/**************************************************/

int getRooInfo(){
  fstream roofile;
  bool closed = true;
  int i, j=0;

  while(closed){
    roofile.open("ini.roo", ios::in);
    if( roofile.is_open() ) closed = false;
    else {
      roofile.close();
      sleep(1);
    }  
  }

  roofile >> start; 
  roofile >> offset; 
  roofile >> distbits; 

  if(rank != 2){
    for(i=0; i<64; i++){
      roofile >> jumpdistance[i]; // = to_ZZ(token);
    }
  }
  else{
    for(i=0; i<64; i++){
      roofile >> jumpdistance[i]; // = to_ZZ(token);
      roofile >> jumpdistance1[i];
    }
  }

  roofile >> gz[0];
  roofile >> gz[1];
  roofile >> gz[2];
  roofile >> gz[3];
  roofile >> gz[4];
  roofile >> gz[5];

  for(i=0; i<64; i++){
    roofile >> jumpsz[6*i];
    roofile >> jumpsz[6*i+1];
    roofile >> jumpsz[6*i+2];
    roofile >> jumpsz[6*i+3];
    roofile >> jumpsz[6*i+4];
    roofile >> jumpsz[6*i+5];
  }

  roofile.close();

  // Create the distinguished point filenames for this kangaroo.
  std::stringstream ss;
  lastdpname = (char *)malloc(17*sizeof(char));

  ss << "last_dp_" << rootype << roonum << ".roo";
  ss >> lastdpname;

  std::stringstream ss2;
  alldpname = (char *)malloc(16*sizeof(char));
  
  ss2 << "all_dp_" << rootype << roonum << ".roo";
  ss2 >> alldpname;
  return(1);
}

/********************************************************/
/* This function reads in the data from the input file. */
/********************************************************/

int getinput(){
  FILE *inputptr;
  char thisLine[128];
  char token[64], value[64];
  int dD, dG, dH, rand;
  int coefg[10];
  int coefh[10];
  int i;
  vec_pair_ZZ_pX_long factors;

  SetSeed( to_ZZ( time(NULL) ) );

  inputptr = fopen("rooinput.txt", "r");

  if(inputptr == NULL){
    cout<<"Error: Could not open rooinput.txt.\n\n";
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
	f = G;
      }
      else if(dH == 1){
	SetX(H);
	f = LeftShift(G, 2);
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

	f = G*sqr(H);
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
      f = G;
    }
    else if(dH == 1) {
      SetX(H);
      
      // Make sure that G and H are relatively prime.
      // And that G is squarefree.
      if(IsZero(coeff(G,0)) || !IsOne(GCD(G, diff(G))) ){
	while(IsZero(coeff(G,0)) || !IsOne(GCD(G, diff(G))) )
	  G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
      }
      f = LeftShift(G, 2);
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

      f = G*sqr(H);
    }
    
  }
  else{
    // Make sure that (G,H) = 1 and G is squarefree.
    if( !IsOne(GCD(G,H)) || !IsOne(GCD(G, diff(G))) ){
      //cout<<"Generating new G(x). Original was either not squarefree or not coprime to H(x).\n\n";
    }
    while( !IsOne(GCD(G,H))  || !IsOne(GCD(G, diff(G))) ){
      G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
    }
    f = G*sqr(H);
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
  if ((dG+2*dH)%3) {
    rank = 0;
    genus = dG + dH -1;
  }
  else {
    genus = dG + dH - 2;
    if(q%3 == 1)
      rank = 2;
    else if(isCube(LeadCoeff(f)))
      rank = 1;
    else
      rank = 0;
  }

  alpha = alphahat();
  if(rank == 1)
    tauratio = tau3();
  if(rank == 2)
    tauratio = tau5();

  fclose(inputptr);

  /********************************/
  /* Get the number of kangaroos. */
  /********************************/

  MPI_Comm_size(MPI_COMM_WORLD, &m);
  m--;
  m-=(m%2); // the number of kangaroos must be even.
  mt = mw = m/2;

  return(1);
 
}

int main(int argc, char *argv[]){
  int found=0;
  int waittime;
  bool closed = true;
  fstream roofile, dpfile, gofile, sfile;
  int i;
  int dptotal = 0, iroodp;
  ZZ gsjumptotal = Zero, bsjumptotal = Zero; 
  ZZ iroogsjump = Zero, iroobsjump = Zero;
  ZZ R, icn;
  RR alphaactual;
  char *filename;
  time_t setk, timek, timeR;
  char *token;
  int rooid;
  int go_signal_out = 1, go_signal_in;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rooid);
  
  /*******************************/
  /* Roo ID 0  = Master Kangaroo */
  /*******************************/

  if(rooid == 0){
    /**************************/
    /* Read in the data file. */
    /**************************/

    if(!getinput()) {
      // If the input could not be obtained, kill the program.
      while(closed){
        roofile.open("solution.txt", ios::out);
        if( roofile.is_open() ) closed = false;
        else  {
          roofile.close();
          sleep(1);
        }  
      }
      roofile<<"Solution not found."<<endl;
      
      roofile.close(); 
      MPI_Finalize();
      
    }
      
    /************************************************************/
    /* Step 1: Compute an approximation E of the class number h */
    /*         and an integer L such that |h - E| < L^2.        */
    /************************************************************/

    approxh();

    /*********************************************/
    /* Step 2: Find elements of small order      */
    /*         then get the initialization info. */
    /*********************************************/
    smallOrder0(a, l);

    if(rank == 0){      
      getRooJumps0();
    }
    else if (rank == 1){
      init(2*genus);
      getRooJumps1();
    }
    else {
      init(2*genus);
      getRooJumps2();
    }

    /***********************************/
    /* Step 3: Send out the kangaroos. */
    /***********************************/

    for(i=1; i<=m; i++)
      MPI_Send(&go_signal_out, 1, MPI_INT, i, 1, MPI_COMM_WORLD);

    // Compute the delay between checking the distinguished point file.
    // This will be calculated so that we expect roughly 10 checks on 
    // average, or once a day, whichever is less.
    // This assumes 1000 jumps are performed each second.
    // Better timing data should be taken to check this.long rate;
    long rate;
    if( rank == 0 ){
      if(genus == 3) rate = 1500;
      else if(genus == 4) rate = 1100;
      else rate = 1000;
      
      waittime = (int)to_long(4*L/RoundToZZ(to_RR(m*40*rate)*sqrt(to_RR(l))));
    }
    else{
      if(genus == 3) rate = 570;
      else if(genus == 4) rate = 300;
      else rate = 300;
      if(rank == 1)
	waittime = (int)to_long(4*L/RoundToZZ(to_RR(m*40*rate)*sqrt(to_RR(2*tauratio-1))));
      else
	waittime = (int)to_long(4*L/RoundToZZ(to_RR(m*40*rate)*sqrt(to_RR(l)*to_RR(2*tauratio-1))));
    }
    
    if(waittime > 3600)//86400)
      waittime = 3600;//86400;
    if(waittime < 2)
      waittime = 2;
    
    // Set the expected size of the dp array.
    size = to_long(RightShift(to_ZZ(waittime)*to_ZZ(m)*10*rate, distbits));
    dp.SetLength(size);
    dists.SetLength(size);
    
    setk = time(NULL);

    /****************************************************************/
    /* Read in the distinguished point file every waittime seconds. */
    /****************************************************************/
    while(!found){
      sleep(waittime);
      found = readDPfile();
    }
    timek = time(NULL);

    // In the unit rank 1 and 2 cases, extract the regulator.
    if(rank == 1){
      
      R = extract(order, One, factors, exponents, 0);
      timeR = time(NULL);
    }
    else if (rank == 2){
      R = extract2(order, One, factors, exponents, 0);
      timeR = time(NULL);
      //if(R != order){
      icn = order/R;
      //latticeRoo(icn, R, e10, e12, e20, e22);
      //}
    }

    /************************************************************/
    /* A solution has been found! Write the solution to a file. */
    /************************************************************/

    while(closed){
      roofile.open("solution.txt", ios::out);
      if( roofile.is_open() ) closed = false;
      else  {
        roofile.close();
        sleep(1);
      }  
    }
    
    roofile<<"G = "<<G<<endl;
    roofile<<"H = "<<H<<endl;
    
    roofile<<"The divisor class number of F_"<<q<<"(C), C: Y^3 = G*H^2 = "<<f<<" is:"<<endl<<endl;
    roofile<<"h   = "<<order;
    if(rank == 0)
      roofile<<"\n\n";
    else{
      roofile<<" = ";
      for(i=0; i < factors.length()-1; i++){
	if(exponents[i] == 1)
	  roofile<<factors[i]<<" * ";
	else
	  roofile<<factors[i]<<"^"<<exponents[i]<<" * ";
      }
      if(exponents[i] == 1)
	roofile<<factors[i];
      else
	roofile<<factors[i]<<"^"<<exponents[i];
      if(!ProbPrime(factors[i]))
	roofile<<" - Last factor not prime"<<endl;
      else
	roofile<<endl;
      if(rank == 1){
	roofile<<"R^S = "<<R<<endl;
	roofile<<"h_x = "<<order/R<<"\n\n";
      }
      else { 
	if(order == R){
	  roofile<<"R   = "<<R<<endl;
	  roofile<<"h_x = "<<icn<<"\n\n";
	}
	else{
	  bool flag = true;
	  for(i=0; i < factors.length()-1; i++){
	    if((exponents[i] > 1) && (IsZero(icn%factors[i]))){
	      flag = false;
	    }
	  }
	  if(flag){
	    R*=icn;  icn = One;
	    roofile<<"R  = "<<R<<endl;
	    roofile<<"h_x = "<<icn<<"\n\n";
	  }
	  else{
	    roofile<<R<<" | R"<<endl;
	    roofile<<"h_x | "<<icn<<"\n\n";
	  }	
	}
      }
    }
    roofile<<"Phase 1 time = "<<phase1time<<" seconds "<<endl;
    roofile<<"Jump time    = "<<timek-setk<<" seconds\n"; 
    if(rank >= 0)
      roofile<<"Phase 4 time = "<<timeR-timek<<" seconds"<<endl;
    roofile<<endl;
    roofile<<"Approximation E = "<<E<<endl;
    roofile<<"Error bound U   = "<<U<<endl;
    alphaactual = to_RR(abs(order-E))/to_RR(U);
    roofile<<"alpha = |h-E|/U = "<<alphaactual<<"\n\n"; 
    
    // Round up the individual distinguished point files.
    token = (char *)malloc(sizeof(char)*100);
    for(i=0; i<mt; i++){
      std::stringstream ss;
      filename = (char *)malloc(sizeof(char)*18);
      ss << "last_dp_t" << i << ".roo";
      ss >> filename;
      // Open the tame distinguished point files.
      dpfile.open(filename,ios::in);
      if( dpfile.is_open() )  {
        dpfile >> token; iroodp = atoi(token); 
        dpfile >> token; iroogsjump = to_ZZ(token);
	if(rank > 0){
	  dpfile >> token; iroobsjump = to_ZZ(token);
        dpfile.close();
        
        dptotal += iroodp;
        gsjumptotal += iroogsjump;
	if(rank > 0)
	  bsjumptotal += iroobsjump;
      }
      else dpfile.close();
    }
    for(i=0; i<mw; i++){ 
      std::stringstream ss;
      filename = (char *)malloc(sizeof(char)*18);
      ss << "last_dp_w" << i << ".roo";
      ss >> filename;
      // Open the wild distinguished point files.
      dpfile.open(filename,ios::in);
      if( dpfile.is_open() )  {
        dpfile>>token; iroodp = atoi(token); 
        dpfile>>token; iroogsjump = to_ZZ(token); 
        if(rank > 0){
	  dpfile >> token; iroobsjump = to_ZZ(token);

        dpfile.close();
        
        dptotal += iroodp;
        gsjumptotal += iroogsjump;
	if(rank > 0)
	  bsjumptotal += iroobsjump;
      }
      else dpfile.close();
    }
    roofile<<"Distinguished points found:      "<<dptotal<<endl;
    roofile<<"Total kangaroo giant step jumps: "<<gsjumptotal<<endl<<endl;
    if(rank > 0){
      roofile<<"Total kangaroo baby step jumps:  "<<bsjumptotal<<endl;
      roofile<<"Total kangaroo jumps:            "<<gsjumptotal+bsjumptotal<<endl;
    }
    roofile<<endl;
    if(rank == 0){
      ZZ expgsjumps = 4*RoundToZZ(to_RR(L)/sqrt(to_RR(l))) + m*power2_ZZ((long)distbits);
      roofile<<"Expected Distinguished Points:   "<<RightShift(expgsjumps,distbits)<<endl;
      roofile<<"Expected total jumps:            "<<expgsjumps<<endl;
    }
    else if (rank == 1){
      RR expgsjumps = 4.0*(to_RR(L)/sqrt(to_RR(2.0*tauratio-1.0)));// + to_RR(m*power2_ZZ((long)distbits))/tauratio;
      ZZ avgjump = CeilToZZ((double)m*sqrt(alpha*to_RR(U)*(2.0*tauratio-1.0)) - 2.0*(tauratio-1.0));
      roofile<<"Expected Distinguished Points:   "<<RightShift(RoundToZZ(expgsjumps),distbits)<<endl;
      roofile<<"Expected total giant step jumps: "<<RoundToZZ(expgsjumps)<<endl;
      roofile<<"Expected total baby step jumps:  "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
      roofile<<"Expected total jumps:            "<<RoundToZZ(tauratio*expgsjumps)<<endl<<endl;
      
      roofile<<"For this example:"<<endl;
      
      expgsjumps = 2.0*((double)m*alphaactual*to_RR(U)/(to_RR(avgjump)+2.0*(tauratio-1.0)) + to_RR(avgjump)/((double)m*(2.0*tauratio-1.0)));// + m*to_RR(power2_ZZ(distbits))/tauratio;
      
      roofile<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
      roofile<<"Expected total giant step jumps: "<<RoundToZZ(expgsjumps)<<endl;
      roofile<<"Expected total baby step jumps:  "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
      roofile<<"Expected total jumps:            "<<RoundToZZ(tauratio*expgsjumps)<<endl;
    }
    else{
      RR expgsjumps = 4.0*(sqrt(alpha*to_RR(U)/(to_RR(l)*(2.0*tauratio-1.0)))) + to_RR(m*power2_ZZ((long)distbits));
      ZZ avgjump = CeilToZZ((double)m*sqrt(alpha*to_RR(U*l)*(2.0*tauratio-1.0))/2.0 - to_RR(l)*(tauratio - 1.0));
      roofile<<"Expected Distinguished Points:   "<<RightShift(RoundToZZ(expgsjumps),distbits)<<endl;
      roofile<<"Expected total giant step jumps: "<<RoundToZZ(expgsjumps)<<endl;
      roofile<<"Expected total baby step jumps:  "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
      roofile<<"Expected total jumps:            "<<RoundToZZ(tauratio*expgsjumps)<<endl<<endl;
      
      roofile<<"For this example:"<<endl;
      
      expgsjumps = ((double)m*alphaactual*to_RR(U)/(to_RR(avgjump)+to_RR(l)*(tauratio-1.0)) + 4.0*to_RR(avgjump)/((double)m*to_RR(l)*(2.0*tauratio-1.0))) + m*to_RR(power2_ZZ(distbits));
      
      roofile<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
      roofile<<"Expected total giant step jumps: "<<RoundToZZ(expgsjumps)<<endl;
      roofile<<"Expected total baby step jumps:  "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
      roofile<<"Expected total jumps:            "<<RoundToZZ(tauratio*expgsjumps)<<endl;
    }
    roofile<<endl;
    roofile<<"Wait time between solution checks: "<<waittime<<" seconds\n\n";
    roofile<<"Used "<<m<<" kangaroos.\n";
    roofile<<"There were "<<mt - m/2<<" collisions among the tame kangaroos.\n";
    roofile<<"There were "<<mw - m/2<<" collisions among the wild kangaroos.\n";
        
    roofile.close();
    
    /*****************************/
    /* Terminate all other work. */
    /*****************************/

    MPI_Finalize();
  }

  /**************************/
  /* Roo ID != 0 = Kangaroo */
  /**************************/

  else {
    MPI_Status rc;

    /****************************/
    /* Initialize the kangaroo. */
    /****************************/

    if(rooid%2)
      rootype = 't';
    else 
      rootype = 'w';
    
    roonum = (rooid-1-((rooid-1)%2))/2;
    
    // Wait for the getS.roo file before doing Phase 1 work.
    MPI_Recv(&go_signal_in, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &rc);

    while(closed){
      sfile.open("getS.roo", ios::in);
      if(sfile.is_open()) closed = false;
      else{
        sfile.close();
        sleep(1);
      }
    }
    sfile>>q;
    ZZ_p::init(q);
    sfile>>G;
    sfile>>H;
    sfile>>m;
    sfile>>lambda;
    sfile.close();
    
    f = G*sqr(H);
    
    if((roonum < 0) || (roonum >= m/2)){
      cout<<"Bad kangaroo number 0 <= i < m/2.\n\n";
      return(0);
    }
    
    closed = true;
    
    // Get the prime splitting info for this block.
    getSplitting();

    /****************************************/
    /* Wait for the go signal before going. */
    /****************************************/

    MPI_Recv(&go_signal_in, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &rc);

    /*********************************/
    /* Get information from roo.ini. */
    /*********************************/
    
    getRooInfo();
    if(rank == 0)
      kangaroo0();
    else if(rank == 1){
      init(2*genus);
      kangaroo1();
    }
    else{
      init(2*genus);
      smallOrder0(a, l);
      kangaroo2();
    }
    MPI_Finalize();
  }
}
