/* 
   classnumber.cc

   This file contains functions to compute the divisor class number, h, 
   of any purely cubic function field. It also contains routines to       
   determine the regulator, R, of a cubic function field of positive unit 
   rank via computing in the infrastructure. Included are routines for
   applying the Baby Step-Giant Step algorithm of Shanks and also the 
   Kangaroo method of Pollard.
   For the functions, we will use results and algorithms from:

   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)
   
   
   How to create the executable:
   --------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download cubic_ideal.cpp, cubic_ideal.h and makefile into 
      ../CubicIdeal.
   Compile libcubic via make in ../CubicIdeal.
   Download infrastructure_ideal.cc, infrastructure_ideal.h and makefile 
      into ../Infideal.
   Compile libinfideal via make in ../Infideal.
   Download hbounds.cc, hbounds.h and makefile 
      into ../Hbounds.
   Compile libhbounds via make in ../Hbounds.
   Download classnumber.cc, classnumber.h and makecn 
      into ../PointCount.
   Compile classnumber via ./makecn in ../PointCount.
   
   Make sure that the paths to the NTL libraries are
   specified correctly in the makefiles.
   
   Known Problems: 
   ---------------
   For unit rank 2 - Test Kangaroo method in Phases 3 and 4.
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   4. Organize this file better and possibly even split it up into several
      files that are more easily managed.

   Author: 
   Eric Landquist

*/

#include "classnumber.h"

NTL_CLIENT

// Baby-step, Giant-step algorithm for unit rank 0.
// Input: a, stepLen - h = a mod stepLen.
// Output: Prints h, |h-E|, and |h-E|/U, along with timing data.
// See Alg. 6.2.4 of [L09].
void bsgs0(ZZ &a, ZZ& stepLen){
  int found = 0, genTrials=0, test2=1; // test1=200;
  long hashsize, place1, place2;
  cubic_ideal g; // A subgroup generator.
  cubic_ideal z, temp, last, zinv;
  cubic_ideal gs; // The giant step.
  cubic_ideal I;  // An identity element.
  ZZ start, end = E+U; // = E-L2, end=E+L2;
  ZZ order, order2, giantstep;
  ZZ ordg;
  ZZ numSteps; // = CeilToZZ(to_RR(L)/sqrt(to_RR(factor)));
  ZZ k = ZZ(), percent, r;
  time_t set, bstime, gstime;
  int pComplete = 0;
  double n = tau1(); // n = #time for one giant step/#time for one baby step.

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

  /**************/
  /* BABY STEPS */
  /**************/
  order0=order1=ZZ();
  cout<<"Phase 3: Making "<<numSteps<<" Baby Steps ... \n";
  set = time(NULL);
  do {
    // Initialize the baby step hash table.

    for(long i=0; i<hashsize; i++){
      baby_loc[i] = N1;
      clear(baby1[i]);
      clear(baby2[i]);
    }

    // Get a random (sub)group generator.
    random(g);

    ordg = makeBabySteps(g, baby1, baby2, baby_loc, hashsize, numSteps, a, stepLen, 1);

    // Check other generators to see if this is the order.
    // Or at least try to find factors of the order.
    if(!IsZero(ordg)){
      genTrials = 0;
      while(genTrials < test2){
	random(z);
	temp = z^ordg;
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
      
	  ordg = abs(order0-order1);
	  temp = g^ordg;
	  // We found a factor.
	  if(IsOne(temp.s) && !IsZero(ordg)){
	    a = ZZ();
	    stepLen = ordg;
	    numSteps = CeilToZZ(to_RR(L)*sqrt(n)/sqrt(to_RR(stepLen)));
	    numSteps += (numSteps%2);
	  }
	  else
	    ordg = stepLen;
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
  giantstep = stepLen*numSteps;
  gs = g^giantstep;
  z = cubic_ideal(gs);

  percent = U/(100*giantstep);
  if(percent==0) percent=1;
  pComplete = 0;

  do{
    do{
      place1 = search(z, baby1, baby2, baby_loc, hashsize);

      if(place1 >= 0)
	found = 1;
      
      // Search the hash table for the inverse of z.
      inverse(zinv, z);
      place2 = search(zinv, baby1, baby2, baby_loc, hashsize);
      
      if(place2 >= 0)
	found = 2;
      
      if(IsZero(k%percent)){
	printf("\r");
	printf("Giant steps are %d %% complete. ", pComplete);

	cout<<flush;
	pComplete++;
	if( pComplete == 100){
	  place2=1;
	}
      }

      temp = cubic_ideal(z);
      z = temp*gs;
      k++;
      if(k==numSteps){
	cout<<"Went overboard!\n\n";
	return;
      }
    }while((place1 < 0) && (place2 < 0));    

    gstime = time(NULL);
    start = E + a - (E%stepLen);

    if(place1 >= 0){
      order = start - giantstep*k + stepLen*baby_loc[place1];
    }
    else{
      order = start + giantstep*k + stepLen*baby_loc[place2];
    }

    cout<<"Match found.\n"<<"Order = "<<order<<endl;

    // Test several elements to see if this is the order.
    genTrials = 0;
    while(genTrials < test2){
      random(g);
      temp = g^order;
      if(temp != I){
	cout<<"Found the wrong order.\n\n";
	found = 0;
	break;
      }
      else
	genTrials++;
    }
    if((place1 >= 0) && (place2 >= 0) && (!found)){
      found = 2;
      order = start + giantstep*k + baby_loc[place2];
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
  } while(!found);

  cout<<endl;
  cout<<k<<" Giant steps completed in "<<gstime-bstime<<" seconds.\n\n";

 end:

  cout<<"The order is h = "<<order<<".\n\n";
  //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<"  |h-E| = "<<abs(order-E)<<endl; 
  //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  //cout<<"|h-E2|/L3^2 = "<<to_RR(abs(order-E))/to_RR(L32)<<endl;
  cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<endl;

  return;
}

// Function to make the baby step set for the unit rank 0 case.
// Input: g: A (sub)group generator.
//        baby1,2: A hash table to store baby steps.
//        loc: An array that gives the power of g stored in the array.
//        hashsize: the size of the hashtable.
//        steps: number of baby steps to compute.
//        a: the equivalence class of the baby steps.
//        stepLen: The length of a baby step.
// Return: either the order of g, or 0 if the identity was not encountered. 
// See Alg. 6.2.9 of [L09].
ZZ makeBabySteps(cubic_ideal &g, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize, ZZ& steps, ZZ &a, ZZ& stepLen, int verbose){
  // The identity and a random (sub)group generator.
  cubic_ideal h, z, step;
  ZZ i = to_ZZ(1), start; 
  ZZ steps0;
  ZZ percent = steps/100;
  if(percent==0)
    percent++;
  int pComplete = 0;

  // We'll be working in equivalence classes a mod stepLen.
  start = E + a - (E%stepLen);
  h = g^start;
  step = g^stepLen;

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
      order0 = start+i*stepLen;
      
      return (order0);
    }
    
    // Otherwise, we insert the value into the hash table.
    insert(h, baby1, baby2, loc, i, hashsize); 
    // Make the next baby step.
    
    z = cubic_ideal(h);
    h = z*step;
  }

  return ZZ();
}

// The Kangaroo Method for unit rank 0.
// Input: a, stepLen - h = a mod stepLen.
// Output: Prints h, |h-E|, and |h-E|/U, along with timing data.
// See Alg. 6.2.9 of [L09].
void kangaroo0(ZZ &a, ZZ& factor){
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
  cubic_ideal jumps[64];
  ZZ jumpdistance[64];
  ZZ jumpsum = ZZ();
  
  ZZ tamedist, wilddist = ZZ(), matchdist = ZZ();

  ZZ jumptotal = ZZ(); // The number of jumps for each kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  long hashsize;
  ZZ order; // The order of the group.
  
  // Set the average jump distance, the jump distances, 
  // and the starting position for tame.
  RoundToZZ(avgjump, to_RR(L)*sqrt(to_RR(factor)));
  
  avgjump-=(avgjump%factor);
  randjump = 2*avgjump/factor;
  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = factor*(RandomBnd(randjump)+1);
      
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*avgjump) || (jumpsum > (i+1)*avgjump) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = factor*(RandomBnd(randjump)+1);
      jumpsum += jumpdistance[i];
    }
    
    jumpdistance[63] -= jumpdistance[i];
  }
  jumpdistance[63] += 64*avgjump;

  
  start = E - E%factor + a;

  distbits = (NumBits(avgjump)-1)/2;
  distbits += (distbits%2);
  distbits2 = distbits/2;

  // Set the hashtables and their size.
  hashsize = NextPrime(to_long(10*RightShift(avgjump, distbits)));
  hashbits = NumBits(hashsize);

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

  for(i=0; i<64; i++){
    jumps[i] = g^jumpdistance[i];
  }

  // Set the tame kangaroo. (The wild kangaroo is initialized at I.)
  tame = g^start;
  tamedist = start;
  //cout<<"BOING! BOING! BOING! BOING! BOING! BOING! BOING! BOING! BOING! \n\n";

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
  //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<" |h-E|  = "<<abs(order-E)<<endl; 
  //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  //cout<<"|h-E2|/L3^2 = "<<to_RR(abs(order-E))/to_RR(U)<<endl;
  cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<endl;
  //if(IsOne(factor)){
  //  cout<<"Expected kangaroo jumps:       "<<4*L<<endl;
  //  cout<<"Expected distinguished points: "<<RightShift(L, distbits-2)<<endl;
  // }
  //else{
  cout<<"Expected kangaroo jumps:       "<<4*avgjump/factor<<endl;
  cout<<"Expected distinguished points: "<<RightShift(avgjump/factor, distbits-2)<<endl;
  //}
  cout<<"Total kangaroo jumps:          "<<jumptotal<<endl;
  cout<<"Total distinguished points:    "<<distpoints<<endl;
}

// Kangaroo method for unit rank 1.
// Input: The function field.
// Output: Prints h, |h-E|, and |h-E|/U, along with timing data.
// See Alg. 6.2.11 of [L09].
void kangaroo1(){
  int i, found = 0;
  int distbits; // The number of zero bits to determine 
                // if an ideal is distinguished.
  int distbits2;
  int dbadd = 4;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  ZZ avgjump;  // The average of the jumps.
  ZZ randjump; // The upper bound on the random jump distance.
  infrastructure_ideal tame, wild;
  infrastructure_ideal temp;
  infrastructure_ideal jumps[64];
  ZZ jumpdistance[64];
  ZZ jumpsum = ZZ();
  
  long headwind;

  ZZ matchdist = ZZ();

  ZZ bsjumptotal = ZZ(); // The number of baby step  jumps for each kangaroo.
  ZZ gsjumptotal = ZZ(); // The number of giant step jumps for each kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  long hashsize;
  ZZ order; // The order of the group.
  ZZ tempnum;
  double tauratio = tau3(); // T_G/T_B - giant step to baby step ratio.
  time_t set, rootime, checktime;
  ZZ Stau;

  cout<<"Initializing ...";
  cout<<flush;

  // Initialize the infrastructure_ideals.
  reduce_basis(tame);
  reduce_basis(wild);
  
  Stau = RoundToZZ(to_RR(q)/tauratio);

  if(q < 1200)
    dbadd-=2;
  if(q < 120)
    dbadd-=2;

  // Set the average jump distance, the jump distances, 
  // and the starting position for tame.
  avgjump = RoundToZZ(2.0*sqrt(alpha*to_RR(U)*(2.0*tauratio-1.0)) - 2.0*(tauratio-1.0));

  // Get the headwind.
  headwind1(headwind);

  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 2 + headwind);
    
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*(avgjump)) || (jumpsum > (i+1)*(avgjump)) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 2 + headwind);
      jumpsum += jumpdistance[i];
    }
    
    jumpdistance[63] -= jumpdistance[i];
  }
  jumpdistance[63] += (64*(avgjump) + headwind + 32);

  // Set the number of bits for a distinguished ideal.
  distbits = (NumBits(avgjump)-1)/2 - 2 + dbadd;

  distbits += (distbits%2);
  distbits2 = distbits/2;

  // Set the hashtables and their size.
  hashsize = NextPrime(to_long(10*RightShift(RoundToZZ(tauratio*to_RR(L)), distbits-2)));
  hashbits = NumBits(hashsize);

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
  
  // Set the tame kangaroo. (The wild kangaroo is initialized at I.)
  tempnum = 2*E;
  below(tempnum, tame);
  
  for(i=0; i<64; i++){
    below(jumpdistance[i], jumps[i]);
  }

  RR expgsjumps = 4.0*sqrt(alpha*to_RR(U)/(2.0*tauratio-1.0)) + to_RR(2*power2_ZZ((long)distbits))/tauratio;

  cout<<" Done. Beginning the jumping."<<endl;
  cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  //cout<<"Expected equiv. GS: "<<to_RR(RoundToZZ(2.0*(      2.0*alpha*to_RR(U)/to_RR(avgjump) + to_RR(avgjump)/(2.0*(2.0*to_RR(tau)-1.0))) + to_RR(power2_ZZ(distbits+1))/((double)tau) )) + to_RR(RoundToZZ(2.0*(tau-1)*(2.0*alpha*to_RR(U)/to_RR(avgjump) + to_RR(avgjump)/(2.0*(2.0*to_RR(tau)-1.0))) + to_RR(power2_ZZ(distbits+1))*(((double)tau-1.0)/(double)tau)))/n<<endl;
  //  cout<<"Expected time: "<<to_RR(to_RR(RoundToZZ(2.0*(      2.0*alpha*to_RR(U)/to_RR(avgjump) + to_RR(avgjump)/(2.0*(2.0*to_RR(tau)-1.0))) + to_RR(power2_ZZ(distbits+1))/((double)tau) )) + to_RR(RoundToZZ(2.0*(tau-1)*(2.0*alpha*to_RR(U)/to_RR(avgjump) + to_RR(avgjump)/(2.0*(2.0*to_RR(tau)-1.0))) + to_RR(power2_ZZ(distbits+1))*(((double)tau-1.0)/(double)tau)))/n)/570.0<<" seconds"<<endl;
  //cout<<"BOING! BOING! BOING! BOING! BOING! BOING! BOING! BOING! BOING! \n\n";
  cout<<"Distinguished points expected:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
  set = time(NULL);
  //cout<<"Initializations completed in "<<set-timingtest<<" seconds"<<endl;
  while(!found){
    // FOR THE TAME KANGAROO:
    
    /************************************************/
    /* Make baby steps until L(tame)(0) < [q/tau_3]. */
    /************************************************/
    while(rep(coeff(tame.d, 0)) >= Stau){
      baby_step_r1(tame, temp);
      tame = temp;

      bsjumptotal++;
    }
    
    // I'll only check for DPs in the set S_tau.
    /***********************************/
    /* Check for distinguished points. */
    /***********************************/
    roocoeff = rep(coeff(tame.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(tame.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(tame, 0, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, matchdist, hashsize); 
	distpoints++;
        
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(tame.d0-matchdist)/2;
	  if(check_inf(order))
	    goto finish;
	  else
	    found = 0;
	}
      }
    }
      
    // FOR THE WILD KANGAROO:
    /************************************************/
    /* Make baby steps until L(wild)(0) < [q/tau_3]. */
    /************************************************/
    while(rep(coeff(wild.d, 0)) >= Stau){
      baby_step_r1(wild, temp);
      wild = temp;

      bsjumptotal++;
    }
    
    /***********************************/
    /* Check for distinguished points. */
    /***********************************/
    roocoeff = rep(coeff(wild.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(wild.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(wild, 1, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, matchdist, hashsize); 
	distpoints++;
	
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(wild.d0-matchdist)/2;
	  if(check_inf(order))
	    goto finish;
	  else
	    found = 0;
	}
      }
    }
    
    // Make the jumps.
    i = vmap(tame);
    temp = tame;
    tame = temp*jumps[i];

    i = vmap(wild);
    temp = wild;
    wild = temp*jumps[i];

    gsjumptotal+=2;     
  }

 finish:

  rootime = time(NULL);

  cout<<"\n"<<"A multiple of R^S is h_0 = "<<order<<".\n\n"; 
  cout<<"Phase 3: Completed in "<<rootime-set<<" seconds."<<endl;
  
  //cout<<" |h-E|  = "<<abs(order-E)<<endl; 
  //cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<"\n"<<endl;
  
  // Extract R_S as a factor of h_0.
  vec_ZZ factors;
  vec_long exponents;
  
  ZZ R = extract(order, One, factors, exponents, manext);
  checktime = time(NULL);

  cout<<endl;
  cout<<"h_0     = "<<order<<" = ";
  for(i=0; i < factors.length()-1; i++){
    if(exponents[i] == 1)
      cout<<factors[i]<<" * ";
    else
      cout<<factors[i]<<"^"<<exponents[i]<<" * ";
  }
  if(exponents[i] == 1)
    cout<<factors[i];
  else
    cout<<factors[i]<<"^"<<exponents[i];
  if(!ProbPrime(factors[i]))
    cout<<" - Last factor not prime"<<endl;
  else
    cout<<endl;

  cout<<"R^S     = "<<R<<endl;
  cout<<"h^*     = "<<order/R<<endl;
  cout<<endl;

  cout<<"Phase 4: Completed in "<<checktime-rootime<<" seconds."<<endl;

  cout<<"\n";
  RR alphaactual = to_RR(abs(order-E))/to_RR(U);
  //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<"|h-E|   = "<<abs(order-E)<<endl; 
  //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  cout<<"|h-E|/U = "<<alphaactual<<endl;
  cout<<endl;
  //cout<<tau<<" "<<alphaactual<<" "<<2.0*alphaactual*to_RR(U)/to_RR(avgjump)<<" "<<to_RR(avgjump)/(2.0*(2.0*to_RR(tau)-1.0))<<" "<<to_RR(power2_ZZ(distbits))<<endl;
  
  expgsjumps = 2.0*(2.0*alphaactual*to_RR(U)/to_RR(avgjump) + to_RR(avgjump)/(2.0*(2.0*tauratio-1.0))) + to_RR(power2_ZZ(distbits+1))/tauratio;
  
  cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  cout<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
  cout<<"Total baby step kangaroo jumps:  "<<bsjumptotal<<endl;
  cout<<"Total giant step kangaroo jumps: "<<gsjumptotal<<endl;
  cout<<"Total kangaroo jumps:            "<<bsjumptotal+gsjumptotal<<endl;
  cout<<"Total distinguished points:      "<<distpoints<<endl;
  cout<<endl;
}

// Kangaroo method for unit rank 2.
// Input: The function field.
// Output: Prints h, |h-E|, and |h-E|/U, along with timing data.
// A generalization of Alg. 6.2.11 of [L09].
void kangaroo2(){
  int i, found = 0;
  int distbits; // The number of zero bits to determine 
                // if an ideal is distinguished.
  int distbits2;
  int dbadd = 4;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  ZZ avgjump;  // The average of the jumps.
  ZZ randjump; // The upper bound on the random jump distance.
  infrastructure_ideal tame, wild;
  infrastructure_ideal temp, C, D;
  cubic_ideal jumps[64];
  ZZ jumpdistance[64];
  ZZ jumpdistance1[64];
  ZZ jumpsum = ZZ();

  ZZ matchdist = ZZ();

  ZZ bsjumptotal = ZZ(); // The number of baby step  jumps for each kangaroo.
  ZZ gsjumptotal = ZZ(); // The number of giant step jumps for each kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  long hashsize;
  ZZ order; // The order of the group.
  double tauratio = tau5(); // T_G/T_B - giant step to baby step ratio.
  time_t set, rootime, checktime;
  ZZ Stau;
  long dist0, dist1, dist2;
  ZZ wind0, wind2; 

  headwind2(dist0, dist1, dist2);
  wind0 = to_ZZ(dist0);
  wind2 = to_ZZ(dist2);
  
  cout<<"Initializing ...";
  cout<<flush;

  // Initialize the infrastructure_ideals.
  reduce_basis(tame);
  reduce_basis(wild);
  below(E, Zero, tame);

  Stau = RoundToZZ(to_RR(q)/tauratio);

  if(q < 1200)
    dbadd-=2;
  if(q < 120)
    dbadd-=2;

  // Set the average jump distance and the jump distances, 
  avgjump = RoundToZZ(sqrt(alpha*to_RR(U)*(2.0*tauratio-1.0)) - tauratio + 1.0);
  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 1 + wind0);
    
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*(avgjump)) || (jumpsum > (i+1)*(avgjump)) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 1 + wind0);
      jumpsum += jumpdistance[i];
    }
    jumpdistance[63] -= jumpdistance[i];
  }
  jumpdistance[63] += (64*(avgjump) + wind0);

  // Set the number of bits for a distinguished ideal.
  distbits = (NumBits(avgjump)-1)/2 - 2 + dbadd;

  distbits += (distbits%2);
  distbits2 = distbits/2;

  RR expgsjumps = 4.0*sqrt(alpha*to_RR(U)/(2.0*tauratio-1.0)) + to_RR(2*power2_ZZ((long)distbits));

  // Set the hashtables and their size.
  hashsize = NextPrime(to_long(10*RightShift(RoundToZZ(expgsjumps), distbits)));
  hashbits = NumBits(hashsize);

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
   
  for(i=0; i<64; i++){
    below(jumpdistance[i], wind2, temp);
    while(temp.d2 != wind2){
      jumpdistance[i]--;
      below(jumpdistance[i], wind2, temp);
    } 
    inf_to_cubic(temp, jumps[i]);
    jumpdistance[i] = temp.d0;
    jumpdistance1[i] = temp.d1;
  }

  cout<<" Done. Beginning the jumping."<<endl;
  cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  cout<<"Distinguished points expected:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
  set = time(NULL);
  //cout<<"Initializations completed in "<<set-timingtest<<" seconds"<<endl;
  while(!found){
    // FOR THE TAME KANGAROO:
    /************************************************/
    /* Make baby steps until L(tame)(0) < [q/tau_3]. */
    /************************************************/
    while(rep(coeff(tame.d, 0)) >= Stau){
      baby_step_0_r2(tame, temp);
      tame = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to 0 as possible.
      if(tame.d2 < 0){
	reduce_basis2(tame);
        while(tame.d2 <= 0){
          C=tame;
          baby_step_2_r2(C, tame);
        }
        reduce_basis(C); 
        tame = C; 
        if(tame.d2 < 0){
          do{
            baby_step_0_r2(tame, C); 
            tame = C;
            reduce_basis2(C);
            while(C.d2 <= 0){
              D = C;
              baby_step_2_r2(D, C);
            }
            tame = D;
            reduce_basis(tame);
          } while(tame.d2 < 0); 
        }   
      }
      bsjumptotal++;
    }

    // I'll only check for DPs in the set S_tau.
    // Check if the ideal is distinguished.
    roocoeff = rep(coeff(tame.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(tame.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(tame, 0, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, matchdist, hashsize); 
	distpoints++;
        
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(tame.d0-matchdist);
	  if(check_inf2(order) && checkPrimes(order))
	    goto finish;
	  else
	    found = 0;
	}
      }
    }
        
    // FOR THE WILD KANGAROO:
    // Make baby steps until s(0) = 0 (mod tau).
    /************************************************/
    /* Make baby steps until L(wild)(0) < [q/tau_3]. */
    /************************************************/
    while(rep(coeff(wild.d, 0)) >= Stau){
      baby_step_0_r2(wild, temp);
      wild = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to 0 as possible.
      if(wild.d2 < 0){
	reduce_basis2(wild);
	while(wild.d2 <= 0){
	  C=wild;
	  baby_step_2_r2(C, wild);
	}
	reduce_basis(C); 
	wild = C; 
	if(wild.d2 < 0){
	  do{
	    baby_step_0_r2(wild, C); 
	    wild = C;
	    reduce_basis2(C);
	    while(C.d2 <= 0){
	      D = C;
	      baby_step_2_r2(D, C);
	    }
	    wild = D;
	    reduce_basis(wild);
	  } while(wild.d2 < 0); 
	}   
      }
      bsjumptotal++;
    }
    
    // Check if the ideal is distinguished.
    roocoeff = rep(coeff(wild.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(wild.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(wild, 1, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, matchdist, hashsize); 
	distpoints++;
	
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(wild.d0-matchdist);
	  if(check_inf2(order) && checkPrimes(order))
	    goto finish;
	  else
	    found = 0;
	}
      }
    }
    
    // Make the jumps.
    i = vmap(tame);
    temp = tame;
    giant_step_r2(tame, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    //tame = temp*jumps[i];

    i = vmap(wild);
    temp = wild;
    giant_step_r2(wild, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    //wild = temp*jumps[i];

    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(tame.d2 < 0){
      reduce_basis2(tame);
      while(tame.d2 <= 0){
	C=tame;
	baby_step_2_r2(C, tame);
      }
      reduce_basis(C); 
      tame = C; 
      if(tame.d2 < 0){
	do{
	  baby_step_0_r2(tame, C); 
	  tame = C;
	  reduce_basis2(C);
	  while(C.d2 <= 0){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  tame = D;
	  reduce_basis(tame);
	} while(tame.d2 < 0); 
      }   
    }
    else if(tame.d2 > 0){
      reduce_basis1(tame);
      while(tame.d2 >= 0){
	D=tame;
	baby_step_1_r2(D, tame);
      }
      reduce_basis(D);
      tame = D;
    }
    if(wild.d2 < 0){
      reduce_basis2(wild);
      while(wild.d2 <= 0){
	C=wild;
	baby_step_2_r2(C, wild);
      }
      reduce_basis(C); 
      wild = C; 
      if(wild.d2 < 0){
	do{
	  baby_step_0_r2(wild, C); 
	  wild = C;
	  reduce_basis2(C);
	  while(C.d2 <= 0){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  wild = D;
	  reduce_basis(wild);
	} while(wild.d2 < 0); 
      }   
    }
    else if(wild.d2 > 0){
      reduce_basis1(wild);
      while(wild.d2 >= 0){
	D=wild;
	baby_step_1_r2(D, wild);
      }
      reduce_basis(D);
      wild = D;
    }
    gsjumptotal+=2;     
  }

 finish:

  rootime = time(NULL);

  cout<<"\n"<<"A multiple of R is h = "<<order<<".\n\n"; 
  cout<<"Phase 3: Completed in "<<rootime-set<<" seconds."<<endl;
  
  //cout<<" |h-E|  = "<<abs(order-E)<<endl; 
  //cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<"\n"<<endl;
  
  // Extract R_S as a factor of h_0.
  vec_ZZ factors;
  vec_long exponents;

  ZZ R = extract2(order, One, factors, exponents, manext);
  
  // If R == order, then we know that h_x = 1 and we are done.
  // Otherwise, the regulator is some multiple of R. Determine which here.
  if(R != order){
    ZZ icn = order/R;
    cout<<"1 <= h_x <= "<<icn<<". Searching the infrastructure to determine h_x and R."<<endl;
    latticeSearch(icn, R);
  }
  checktime = time(NULL);

  cout<<endl;
  cout<<"h     = "<<order<<" = ";
  for(i=0; i < factors.length()-1; i++){
    if(exponents[i] == 1)
      cout<<factors[i]<<" * ";
    else
      cout<<factors[i]<<"^"<<exponents[i]<<" * ";
  }
  if(exponents[i] == 1)
    cout<<factors[i];
  else
    cout<<factors[i]<<"^"<<exponents[i];
  if(!ProbPrime(factors[i]))
    cout<<" - Last factor not prime"<<endl;
  else
    cout<<endl;

  cout<<"R     = "<<R<<endl;
  cout<<"h_x   = "<<order/R<<endl;
  cout<<endl;

  cout<<"Phase 4: Completed in "<<checktime-rootime<<" seconds."<<endl;

  cout<<"\n";
  RR alphaactual = to_RR(abs(order-E))/to_RR(U);
  //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<"|h-E|   = "<<abs(order-E)<<endl; 
  //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  cout<<"|h-E|/U = "<<alphaactual<<endl;
  cout<<endl;
  //cout<<tau<<" "<<alphaactual<<" "<<2.0*alphaactual*to_RR(U)/to_RR(avgjump)<<" "<<to_RR(avgjump)/(2.0*(2.0*to_RR(tau)-1.0))<<" "<<to_RR(power2_ZZ(distbits))<<endl;
  
  expgsjumps = (2.0*alphaactual*to_RR(U)/to_RR(avgjump) + 2.0*to_RR(avgjump)/(2.0*tauratio-1.0)) + to_RR(power2_ZZ(distbits));
  cout<<"For this example:"<<endl;
  cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  cout<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  cout<<endl;
  cout<<"Total baby step kangaroo jumps:  "<<bsjumptotal<<endl;
  cout<<"Total giant step kangaroo jumps: "<<gsjumptotal<<endl;
  cout<<"Total kangaroo jumps:            "<<bsjumptotal+gsjumptotal<<endl;
  cout<<"Total distinguished points:      "<<distpoints<<endl;
  cout<<endl;
}

// Kangaroo method for unit rank 2.
// Input: a - an integer
//        l - an integer such that h = a (mod l)
// Output: Prints h, |h-E|, and |h-E|/U, along with timing data.
// A generalization of Alg. 6.2.11 of [L09].
void kangaroo2(ZZ &a, ZZ &l){
  int i, found = 0;
  int distbits; // The number of zero bits to determine 
                // if an ideal is distinguished.
  int distbits2;
  int dbadd = 4;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  ZZ avgjump;  // The average of the jumps.
  ZZ randjump; // The upper bound on the random jump distance.
  infrastructure_ideal temp, tame, wild, C, D;
  cubic_ideal jumps[64];
  ZZ jumpdistance[64];
  ZZ jumpdistance1[64];
  ZZ jumpsum = ZZ();
  
  ZZ matchdist = ZZ();

  ZZ bsjumptotal = ZZ(); // The number of baby step  jumps for each kangaroo.
  ZZ gsjumptotal = ZZ(); // The number of giant step jumps for each kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  long hashsize;
  ZZ order; // The order of the group.
  double tauratio = tau5(); // T_G/T_B - giant step to baby step ratio.
  time_t set, rootime, checktime;
  ZZ Stau;
  long dist0, dist1, dist2;
  ZZ wind0, wind2; 

  if(tauratio <= to_double(l)){
    tauratio = 1.0;
  }
  else{
    tauratio /= to_double(l);
  }

  headwind2(dist0, dist1, dist2);
  wind0 = to_ZZ(dist0);
  wind2 = to_ZZ(dist2);
  
  cout<<"Initializing ...";
  cout<<flush;

  // Initialize the infrastructure_ideals.
  reduce_basis(tame);
  reduce_basis(wild);
  E -= (E%l); E += a;
  below(E, Zero, tame);
  while (tame.d0%l != a){
    E += l;
    below(E, Zero, tame);
  }
  Stau = RoundToZZ(to_RR(q)/tauratio);

  if(q < 1200)
    dbadd-=2;
  if(q < 120)
    dbadd-=2;

  // Set the average jump distance and the jump distances, 
  avgjump = RoundToZZ(sqrt(alpha*to_RR(U*l)*(2.0*tauratio-1.0)) - to_RR(l)*(tauratio - 1.0));
  randjump = (2*avgjump - (genus+1))/l;
  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = l*(RandomBnd(randjump) + genus + 1) + wind0;
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*(avgjump)) || (jumpsum > (i+1)*(avgjump)) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = l*(RandomBnd(randjump) + genus + 1) + wind0;
      jumpsum += jumpdistance[i];
    }
    jumpdistance[63] -= jumpdistance[i];
  }
  jumpdistance[63] += 64*(avgjump);
  jumpdistance[63] -= jumpdistance[63]%l;
  jumpdistance[63] += wind0;


  // Set the number of bits for a distinguished ideal.
  distbits = (NumBits(avgjump)-1)/2 - 2 + dbadd;

  distbits += (distbits%2);
  distbits2 = distbits/2;

  RR expgsjumps = 4.0*sqrt(alpha*to_RR(U)/(to_double(l)*(2.0*tauratio-1))) + to_RR(2*power2_ZZ((long)distbits));

  // Set the hashtables and their size.
  hashsize = NextPrime(to_long(10*RightShift(RoundToZZ(expgsjumps), distbits)));
  hashbits = NumBits(hashsize);

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
   
  // Set the jumps.
  for(i=0; i<64; i++){
    below(jumpdistance[i], wind2, temp);
    while((temp.d0%l != wind0) && (temp.d2 != wind2)){
      jumpdistance[i]+=l;
      below(jumpdistance[i], wind2, temp);
    } 
    inf_to_cubic(temp, jumps[i]);
    jumpdistance[i] = temp.d0;
    jumpdistance1[i] = temp.d1;
  }

  cout<<" Done. Beginning the jumping."<<endl;
  cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  cout<<"Distinguished points expected:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
  set = time(NULL);
  //cout<<"Initializations completed in "<<set-timingtest<<" seconds"<<endl;
  while(!found){
    // FOR THE TAME KANGAROO:
    /************************************************/
    /* Make baby steps until L(tame)(0) < [q/tau_3]. */
    /************************************************/
    while( (rep(coeff(tame.d, 0)) >= Stau) && (tame.d0%l == a) ){
      baby_step_0_r2(tame, temp);
      tame = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to 0 as possible.
      if(tame.d2 < 0){
	reduce_basis2(tame);
	while(tame.d2 <= 0){
	  C=tame;
	  baby_step_2_r2(C, tame);
	}
	reduce_basis(C); 
	tame = C; 
	if(tame.d2 < 0){
	  do{
	    baby_step_0_r2(tame, C); 
	    tame = C;
	    reduce_basis2(C);
	    while(C.d2 <= 0){
	      D = C;
	      baby_step_2_r2(D, C);
	    }
	    tame = D;
	    reduce_basis(tame);
	  } while(tame.d2 < 0); 
	}   
	
      }

      bsjumptotal++;
    }
    
    // I'll only check for DPs in the set S_tau.
    // Check if the ideal is distinguished.
    roocoeff = rep(coeff(tame.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(tame.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(tame, 0, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, matchdist, hashsize); 
	distpoints++;
        
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(tame.d0-matchdist);
	  if(order%l == a){
	    if(check_inf2(order) && checkPrimes(order))
	      goto finish;
	    else
	      found = 0;
	  }
	  else
	    found = 0;
	}
      }
    }
        
    // FOR THE WILD KANGAROO:
    // Make baby steps until s(0) = 0 (mod tau).
    /************************************************/
    /* Make baby steps until L(wild)(0) < [q/tau_3]. */
    /************************************************/
    while( (rep(coeff(wild.d, 0)) >= Stau) && IsZero(wild.d0%l) ){
      baby_step_0_r2(wild, temp);
      wild = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to 0 as possible.
      if(wild.d2 < 0){
	reduce_basis2(wild);
	while(wild.d2 <= 0){
	  C=wild;
	  baby_step_2_r2(C, wild);
	}
	reduce_basis(C); 
	wild = C; 
	if(wild.d2 < 0){
	  do{
	    baby_step_0_r2(wild, C); 
	    wild = C;
	    reduce_basis2(C);
	    while(C.d2 <= 0){
	      D = C;
	      baby_step_2_r2(D, C);
	    }
	    wild = D;
	    reduce_basis(wild);
	  } while(wild.d2 < 0); 
	}   
      }
      bsjumptotal++;
    }
    
    // Check if the ideal is distinguished.
    roocoeff = rep(coeff(wild.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(wild.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(wild, 1, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, matchdist, hashsize); 
	distpoints++;
	
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(wild.d0-matchdist);
	  if(order%l == a){
	    if(check_inf2(order) && checkPrimes(order))
	      goto finish;
	    else
	      found = 0;
	  }
	  else
	    found = 0;
	}
      }
    }
   
    // Make the jumps.
    i = vmap(tame);
    temp = tame;
    giant_step_r2(tame, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    //tame = temp*jumps[i];
    
    i = vmap(wild);
    temp = wild;
    giant_step_r2(wild, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    //wild = temp*jumps[i];

    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(tame.d2 < 0){
      reduce_basis2(tame);
      while(tame.d2 <= 0){
	C=tame;
	baby_step_2_r2(C, tame);
      }
      reduce_basis(C); 
      tame = C; 
      if(tame.d2 < 0){
	do{
	  baby_step_0_r2(tame, C); 
	  tame = C;
	  reduce_basis2(C);
	  while(C.d2 <= 0){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  tame = D;
	  reduce_basis(tame);
	} while(tame.d2 < 0); 
      }   
    }
    else if(tame.d2 > 0){
      reduce_basis1(tame);
      while(tame.d2 >= 0){
	temp=tame;
	baby_step_1_r2(temp, tame);
      }
      reduce_basis(temp);
      tame = temp;
    }
    if(wild.d2 < 0){
      reduce_basis2(wild);
      while(wild.d2 <= 0){
	C=wild;
	baby_step_2_r2(C, wild);
      }
      reduce_basis(C); 
      wild = C; 
      if(wild.d2 < 0){
	do{
	  baby_step_0_r2(wild, C); 
	  wild = C;
	  reduce_basis2(C);
	  while(C.d2 <= 0){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  wild = D;
	  reduce_basis(wild);
	} while(wild.d2 < 0); 
      }   
    }
    else if(wild.d2 > 0){
      reduce_basis1(wild);
      while(wild.d2 >= 0){
	temp=wild;
	baby_step_1_r2(temp, wild);
      }
      reduce_basis(temp);
      wild = temp;
    }
    // Get the jumps back in the correct residue classes.

    while(tame.d0%l != a){
      temp = tame;
      baby_step_0_r2(temp, tame);
    }
    while(!IsZero(wild.d0%l)){
      temp = wild;
      baby_step_0_r2(temp, wild);
    }
    gsjumptotal+=2;     
  }

 finish:

  rootime = time(NULL);

  cout<<"\n"<<"A multiple of R is h = "<<order<<".\n\n"; 
  cout<<"Phase 3: Completed in "<<rootime-set<<" seconds."<<endl;
  
  //cout<<" |h-E|  = "<<abs(order-E)<<endl; 
  //cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<"\n"<<endl;
  
  // Extract R_S as a factor of h_0.
  vec_ZZ factors;
  vec_long exponents;

  ZZ R = extract2(order, One, factors, exponents, manext);
  
  // If R == order, then we know that h_x = 1 and we are done.
  // Otherwise, the regulator is some multiple of R. Determine which here.
  if(R != order){
    ZZ icn = order/R;
    cout<<"1 <= h_x <= "<<icn<<". Searching the infrastructure to determine h_x and R."<<endl;
    latticeSearch(icn, R);
  }
  checktime = time(NULL);

  cout<<endl;
  cout<<"h     = "<<order<<" = ";
  for(i=0; i < factors.length()-1; i++){
    if(exponents[i] == 1)
      cout<<factors[i]<<" * ";
    else
      cout<<factors[i]<<"^"<<exponents[i]<<" * ";
  }
  if(exponents[i] == 1)
    cout<<factors[i];
  else
    cout<<factors[i]<<"^"<<exponents[i];
  if(!ProbPrime(factors[i]))
    cout<<" - Last factor not prime"<<endl;
  else
    cout<<endl;

  cout<<"R     = "<<R<<endl;
  cout<<"h_x   = "<<order/R<<endl;
  cout<<endl;

  cout<<"Phase 4: Completed in "<<checktime-rootime<<" seconds."<<endl;

  cout<<"\n";
  RR alphaactual = to_RR(abs(order-E))/to_RR(U);
  //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<"|h-E|   = "<<abs(order-E)<<endl; 
  //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  cout<<"|h-E|/U = "<<alphaactual<<endl;
  cout<<endl;
  //cout<<tau<<" "<<alphaactual<<" "<<2.0*alphaactual*to_RR(U)/to_RR(avgjump)<<" "<<to_RR(avgjump)/(2.0*(2.0*to_RR(tau)-1.0))<<" "<<to_RR(power2_ZZ(distbits))<<endl;
  
  expgsjumps = (2.0*alphaactual*to_RR(U)/to_RR(avgjump) + 2.0*to_RR(avgjump)/(to_double(l)*(2.0*tauratio-1.0))) + to_RR(power2_ZZ(distbits));
  cout<<"For this example:"<<endl;
  cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ(to_double(l)*(tauratio-1.0)*expgsjumps)<<endl;
  cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  cout<<"Expected kangaroo jumps:         "<<RoundToZZ(to_double(l)*tauratio*expgsjumps)<<endl;
  cout<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  cout<<endl;
  cout<<"Total baby step kangaroo jumps:  "<<bsjumptotal<<endl;
  cout<<"Total giant step kangaroo jumps: "<<gsjumptotal<<endl;
  cout<<"Total kangaroo jumps:            "<<bsjumptotal+gsjumptotal<<endl;
  cout<<"Total distinguished points:      "<<distpoints<<endl;
  cout<<endl;
}

// A hash function for unit rank 0.
// Input: A - a kangaroo  
// Output: An element of {0, 1, ..., 63}

inline int vmap(cubic_ideal &A){
  return(rep(coeff(A.u,0))%64);
}

// A hash function for unit rank 1 and 2.
// Input: A - a kangaroo  
// Output: An element of {0, 1, ..., 63}

inline int vmap(infrastructure_ideal &A){
  return(rep(coeff(A.mu1,0))%64);
}

// Processes distinguished points (traps) for the Kangaroo Method in rank 0.
// Searches the distinguished point arrays for a match, and inserts otherwise.
// Input: roo - the kangaroo
//        torw = 0 if roo is tame
//               1 if roo is wild
//        dpt1,2 - distinguished point arrays for tame roos
//        dptd - distances of the tame distinguished points
//        dpw1,2 - distinguished point arrays for wild roos
//        dpwd - distances of the wild distinguished points
//        distance - distance of roo
//        size - size of dpt1,2 and dpwd1,2.
// Output: match - the distance of the matched roo
// Return: 0 if there is no match 
//         1 if there is.
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

  return 0;

}
// Processes distinguished points (traps) for the Kangaroo Method in ranks 1,2.
// Searches the distinguished point arrays for a match, and inserts otherwise.
// Input: roo - the kangaroo
//        torw = 0 if roo is tame
//               1 if roo is wild
//        dpt1,2 - distinguished point arrays for tame roos
//        dptd - distances of the tame distinguished points
//        dpw1,2 - distinguished point arrays for wild roos
//        dpwd - distances of the wild distinguished points
//        distance - distance of roo
//        size - size of dpt1,2 and dpwd1,2.
// Output: match - the distance of the matched roo
// Return: 0 if there is no match 
//         1 if there is.
int distpoint(infrastructure_ideal &roo, int torw, ZZVec &dpt1, ZZVec &dpt2, ZZVec &dptd, ZZVec &dpw1, ZZVec &dpw2, ZZVec &dpwd, ZZ &match, long size){
  long index;

  // If the roo is tame, search for a match in the wild dp tables.
  if(torw == 0){
    index = roosearch(roo, dpw1, dpw2, dpwd, size);
    // If we didn't find a match, insert the tame kangaroo into the hash table.
    if(index == -1){
      rooinsert(roo, dpt1, dpt2, dptd, size);
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
      rooinsert(roo, dpw1, dpw2, dpwd, size);
      return 0;
    }
    // Sweet! We found a match!
    else{
      match = dptd[index];
      return 1;
    }
  }

  return 0;

}

// The hash function for inserting distinguished points into the hash table
// for unit rank 0, kangaroo method.
// Input: A - the kangaroo
//        size - the size of the hashtable.

inline long roohash(cubic_ideal &A, long size){
  return( ((rep(coeff(A.v,0))<<hashbits) + rep(coeff(A.s,1)))%size );
}

// The hash function for inserting distinguished points into the hash table
// for unit ranks 1 and 2, kangaroo method.
// Input: A - the kangaroo
//        size - the size of the hashtable.
inline long roohash(infrastructure_ideal &A, long size){
  return( ((rep(coeff(A.nu1, 0))<<hashbits) + rep(coeff(A.mu0, 1)))%size );
}

// Function to insert a kangaroo into its hash table, rank 0.
// Input: A - the kangaroo
//        roo1,2 - the hash tables
//        loc - the distance hash table
//        distance - the distance of A
//        size - the size of the hashtable.
inline void rooinsert(cubic_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, ZZ distance, long size){
  long index, j;

  index = roohash(A, size);
      
  // Collision resolution - quadratic probing.
  j = 1;
  while(loc[index] != N1){
    index+=(j*j);
    index%=size;
    j++;
  }
  
  if(j>maxCollision)
    maxCollision = j;
  
  cout<<"Insert at "<<index<<endl;

  // Insert the kangaroo.
  roo1[index] = rep(coeff(A.s,1));
  roo2[index] = rep(coeff(A.s,2));
  loc[index] = distance;
}

// Function to insert a kangaroo into its hash table, ranks 1 and 2.
// Input: A - the kangaroo
//        roo1,2 - the hash tables
//        loc - the distance hash table
//        distance - the distance of A
//        size - the size of the hashtable.
inline void rooinsert(infrastructure_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, long size){
  long index, j;

  index = roohash(A, size);
      
  // Collision resolution - quadratic probing.
  j = 1;
  while(loc[index] != N1){
    index+=(j*j);
    index%=size;
    j++;
  }
  
  if(j>maxCollision)
    maxCollision = j;
  
  // Insert the kangaroo.
  roo1[index] = rep(coeff(A.mu2, 0));
  roo2[index] = rep(coeff(A.nu2, 0));
  loc[index] = A.d0;
}

// Function to search for a kangaroo in the other hash table, rank 0.
// Input: A - the kangaroo
//        roo1,2 - the hash tables
//        loc - the distance hash table
//        hashsize - the size of the hashtable.
// Output: The position of A if it is found.
//         -1 if A is not found. 
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

// Function to search for a kangaroo in the other hash table, ranks 1 and 2.
// Input: A - the kangaroo
//        roo1,2 - the hash tables
//        loc - the distance hash table
//        hashsize - the size of the hashtable.
// Output: The position of A if it is found.
//         -1 if A is not found. 
inline long roosearch(infrastructure_ideal &A, ZZVec &roo1, ZZVec &roo2, ZZVec &loc, long hashsize){
  long j, index;
  int i;

  index = roohash(A, hashsize);

  if(loc[index] == N1)  
    return(-1);
  
  if(roo1[index] == rep(coeff(A.mu2, 0))){
    if(roo2[index] == rep(coeff(A.nu2, 0)))
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
    if(roo1[index] == rep(coeff(A.mu2, 0))){
      if(roo2[index] == rep(coeff(A.nu2, 0)))
	return index;
    }
  }
  return(-1);
}

// The hash function for inserting baby steps into the hash table
// for unit rank 0, BSGS method.
// Input: A - the baby step
//        size - the size of the hashtable.
// Output: The hash value.
inline long hash(cubic_ideal &A, long size){
  long value = 0;
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
    }
  }
  else{ 
    value = ((rep(coeff(A.u, 0)) << 2) + rep(coeff(A.v,0))) % size;
    for(i=1; i<=d; i++){
      value = ((value << (hashbits/2)) + rep(coeff(A.u, i))) % size; 
    }
    return(value);
  }
}

// The hash function for inserting baby steps into the hash table
// for unit ranks 1 and 2, BSGS method.
// Input: A - the baby step
//        size - the size of the hashtable.
// Output: The hash value.
inline long hash(infrastructure_ideal &A, long size){
  long value=0;
  int i, d = deg(A.mu0);
  
  for(i=0; i<=d; i++){
    value = ((value << (hashbits)) + rep(coeff(A.mu0, i))) % size; 
  }
  d = deg(A.mu1);
  for(i=0; i<=d; i++){
    value = ((value << (hashbits)) + rep(coeff(A.mu1, i))) % size; 
  }
  return(value);
}

// Function to insert a baby step into the hash table, rank 0.
// Input: A - the baby step
//        baby1,2 - the hash tables
//        loc - the distance hash table
//        val - the discrete log of A
//        size - the size of the hashtable.
inline void insert(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, ZZ val, long size){
  long index, j;

  // Otherwise, we insert the value into the hash table.
  index = hash(A, size);
  
  // Collision resolution - quadratic probing.
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

// Function to insert a baby step into the hash table, ranks 1 and 2.
// Input: A - the baby step
//        baby1,2 - the hash tables
//        loc - the distance hash table
//        val - the distance of A
//        size - the size of the hashtable.
inline void insert(infrastructure_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long size){
  long index, j;

  // Otherwise, we insert the value into the hash table.
  index = hash(A, size);
  // Collision resolution - quadratic probing.
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
  loc[index] = A.d0;
}

// Function to search for a giant step in the baby step table, rank 0.
// Input: A - the giant step
//        baby1,2 - the baby step hash tables
//        loc - the distance hash table
//        hashsize - the size of the hashtable.
// Output: The position of A if it is found.
//         -1 if A is not found. 
inline long search(cubic_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize){
  long j, index;
  int i;
  ZZ bvalue1, bvalue2;

  index = hash(A, hashsize);
  bvalue1 = baby1val(A);
  bvalue2 = baby2val(A);
  
  if(loc[index] == N1)  
    return(-1);

  if(baby1[index] == bvalue1){ 
    if(IsOne(A.s)){
      if(IsZero(baby2[index]))
	return index;
    }
    else{
      if(baby2[index] == bvalue2)
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
    if(baby1[index] == bvalue1){ 
      if(IsOne(A.s)){
	if(IsZero(baby2[index]))
	  return index;
      }
      else{
	if(baby2[index] == bvalue2) 
	  return index;
      }
    }
  }
  
  return(-1);
}

// Function to search for a giant step in the baby step table, ranks 1 and 2.
// Input: A - the giant step
//        baby1,2 - the baby step hash tables
//        loc - the distance hash table
//        hashsize - the size of the hashtable.
// Output: The position of A if it is found.
//         -1 if A is not found. 
inline long search(infrastructure_ideal &A, ZZVec &baby1, ZZVec &baby2, ZZVec &loc, long hashsize){
  long j, index;
  int i;
  ZZ bvalue1, bvalue2;

  index = hash(A, hashsize);
  bvalue1 = baby1val(A);
  bvalue2 = baby2val(A);
  if(loc[index] == N1)  
    return(-1);

  if(baby1[index] == bvalue1){ 
    if(baby2[index] == bvalue2)
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
    if(baby1[index] == bvalue1){ 
      if(baby2[index] == bvalue2) 
	return index;
    }
  }
  return(-1);
}

// Input: A - a baby or giant step, unit rank 0.
// Output: value: some information from A to store in the hash table.
// For the sake of space, we don't store an entire ideal in memory,
// but only that which would identify it uniquely with high probability.
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
// Input: A - a baby or giant step, unit rank 0.
// Output: value: some information from A to store in the hash table.
// For the sake of space, we don't store an entire ideal in memory,
// but only that which would identify it uniquely with high probability.
inline ZZ baby2val(cubic_ideal &A){

  int i, d = deg(A.s);
  ZZ value = ZZ();
  
  if(IsOne(A.s))
    return(value); 
  else{
    value = rep(coeff(A.s,1));
    for(i=3; i<d; i+=2){
      value = (value << 4) + rep(coeff(A.s, i));
    }
    return value;
  }
}
// Input: A - a baby or giant step, unit ranks 1 and 2.
// Output: value: some information from A to store in the hash table.
// For the sake of space, we don't store an entire ideal in memory,
// but only that which would identify it uniquely with high probability.
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
// Input: A - a baby or giant step, unit ranks 1 and 2.
// Output: value: some information from A to store in the hash table.
// For the sake of space, we don't store an entire ideal in memory,
// but only that which would identify it uniquely with high probability.
inline ZZ baby2val(infrastructure_ideal &A){

  int i, d = deg(A.nu0);
  ZZ value = ZZ();

  if(IsOne(A.nu0))
    return(value); 
  else{
    value = rep(coeff(A.nu0,0));
    for(i=1; i<d; i+=2){
      value = (value << 4) + rep(coeff(A.nu0, i));
    }
    return value;
  }
}
// Function to make the baby step set for the unit rank 1 case.
// Input: baby1,2 - A hash table to store baby steps.
//        baby_deg - An array that stores the distances.
//        hashsize - the size of the hashtable.
// Output: GS - the giant step to use for the giant step portion.
// Return: A multiple of R, if discovered by chance 
//         or 0 (most likely) otherwise. 
// See Alg. 6.2.6 of [L09].
ZZ makeBabySteps1(ZZVec& baby1, ZZVec& baby2, ZZVec& baby_deg, long hashsize, infrastructure_ideal& GS){
  ZZ i;
  int pComplete = 1;
  ZZ percent = L/100;
  ZZ tempnum;
  infrastructure_ideal A, B, C;
  vec_pair_ZZ_pX_long factors;
  
  //E = 0;
  //L = 50;

  tempnum = 2*E;
  below(tempnum, A);
  C = A;

  insert(A, baby1, baby2, baby_deg, hashsize);

  if(IsZero(percent))
    percent++;

  for(i=1; i<L; i++){
    //cout<<A.d<<" "<<A.d0<<endl;
    if(IsZero(i%percent)){
      printf("\r");
      printf("Baby steps are %d %% complete.", pComplete);
      cout<<flush;
      pComplete++;
    }
    baby_step_r1(A, B);
    A = B;

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
  tempnum = A.d0 - C.d0 - 2;
  below(tempnum, GS);

  return ZZ();
}
// Function to make the baby step set for the unit rank 2 case.
// Input: baby1,2 - A hash table to store baby steps.
//        baby_deg - An array that stores the distances.
//        hashsize - the size of the hashtable.
// Output: GS - the giant step to use for the giant step portion.
// Return: A factor of h, if discovered by chance 
//         or 0 (most likely) otherwise. 
// Extension of Alg. 6.2.6 of [L09].
ZZ makeBabySteps2(ZZVec& baby1, ZZVec& baby2, ZZVec& baby_deg, long hashsize, infrastructure_ideal& GS){
  ZZ i;
  int pComplete = 1;
  ZZ percent = L/100;
  infrastructure_ideal A, B, C, D;
  vec_pair_ZZ_pX_long factors;
  ZZ tempnum;
  long dist0, dist1, dist2;
  ZZ wind2; 

  headwind2(dist0, dist1, dist2);
  wind2 = to_ZZ(dist2);

  below(E, Zero, A);
  C = A;
  //cout<<A.d0-E<<" "<<A.d0<<" "<<A.d1<<" "<<A.d2<<" "<<A.d<<endl;
  insert(A, baby1, baby2, baby_deg, hashsize);

  if(IsZero(percent))
    percent++;

  for(i=1; i<L; i++){
    
    if(IsZero(i%percent)){
      printf("\r");
      printf("Baby steps are %d %% complete.", pComplete);
      cout<<flush;
      pComplete++;
    }
    baby_step_0_r2(A, B);
    A = B;
    //if(A.d0-E == 1992){// && (A.d0-E < 2000))
    //  cout<<"1: "<<A.d0<<" "<<A.d1<<" "<<A.d2<<" "<<A.d<<endl;
    //  A.print();
    // }
    insert(A, baby1, baby2, baby_deg, hashsize);
    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(A.d2 < 0){
      while(A.d2 <= 0){
	reduce_basis2(A);
	D=A;	
	baby_step_2_r2(D, A);
	reduce_basis(A);
	//if((A.d0-E > 1985) && (A.d0-E < 2000))
	//  cout<<"2: "<<A.d0<<" "<<A.d1<<" "<<A.d2<<" "<<A.d<<endl;
	insert(A, baby1, baby2, baby_deg, hashsize);
      }
      reduce_basis(D);
      A = D;
      if(A.d2 < 0){
	do{
	  baby_step_0_r2(A, D);
	  //if((D.d0-E > 1985) && (D.d0-E < 2000))
	  //  cout<<"3: "<<D.d0<<" "<<D.d1<<" "<<D.d2<<" "<<D.d<<endl;
	  insert(D, baby1, baby2, baby_deg, hashsize);
	  A = D;
	  while(D.d2 <= 0){
	    reduce_basis2(D);
	    B = D;
	    baby_step_2_r2(B, D);
	    reduce_basis(D);
	    //if((D.d0-E > 1985) && (D.d0-E < 2000))
	    //  cout<<"4: "<<D.d0<<" "<<D.d1<<" "<<D.d2<<" "<<D.d<<endl;
	    insert(D, baby1, baby2, baby_deg, hashsize);
	  }
	  A = B;
	  reduce_basis(A);
	} while(A.d2 < 0);
      }
    }
    
    // In the off-chance that we actually make it around one dimension
    // of the lattice, we've found a factor of h!
    if(IsOne(A.d)){
      GS = infrastructure_ideal(A);
      return(A.d0);
    }
    if(C == A){
      GS = A;
      return(A.d0 - C.d0);
    }
  }
  tempnum = A.d0 - C.d0;
  do {
    below(tempnum, wind2, GS);
    tempnum--;
  } while(GS.d2 != wind2);
  return Zero;
}

// Function to make the baby step set for the unit rank 2 case in Phase 4.
// Input: baby1,2 - A hash table to store baby steps.
//        baby_deg - An array that stores the distances.
//        hashsize - the size of the hashtable.
// Output: GS - the giant step to use for the giant step portion of Phase 4.
void latticeBS(ZZ &bsNum, ZZVec &baby1, ZZVec &baby2, ZZVec &baby_deg, long hashsize, infrastructure_ideal &GS){
  infrastructure_ideal A, B, D;
  ZZ i, percent = bsNum/100;
  int pComplete = 1;
  ZZ tempnum, wind;
  long wind0, wind1, wind2;

  headwind2(wind0, wind1, wind2);
  wind = to_ZZ(wind2);
  
  reduce_basis(A);

  insert(A, baby1, baby2, baby_deg, hashsize);

  if(IsZero(percent))
    percent++;

  for(i=1; i<bsNum; i++){
    if(IsZero(i%percent)){
      printf("\r");
      printf("Baby steps are %d %% complete.", pComplete);
      cout<<flush;
      pComplete++;
    }
    //cout<<i<<"-> "<<A.d0<<" "<<A.d1<<" "<<A.d2<<" "<<A.d<<endl;
    baby_step_0_r2(A, B);
    A = B;
    //if(A.d0 == 5889)
    //  cout<<i<<"-> "<<A.d0<<" "<<A.d1<<" "<<A.d2<<" "<<A.d<<endl;
    
    insert(A, baby1, baby2, baby_deg, hashsize);
    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(A.d2 < 0){
      //if(i==1)
      //cout<<"1: "<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
      while(A.d2 <= 0){
	reduce_basis2(A);
	D=A;	
	baby_step_2_r2(D, A);
	reduce_basis(A);
	//if(A.d0==5889)
	// cout<<"1a: "<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
	insert(A, baby1, baby2, baby_deg, hashsize);
      }
      reduce_basis(D);
      A = D;
      //if(A.d0==5889)
      //cout<<"2: "<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
      if(A.d2 < 0){
	do{
	  baby_step_0_r2(A, D);
	  //if(D.d0==5889)
	  //  cout<<"2a: "<<D.d0<<" "<<D.d1<<" "<<D.d2<<endl;
	  insert(D, baby1, baby2, baby_deg, hashsize);
	  A = D;
	  while(D.d2 <= 0){
	    reduce_basis2(D);
	    B = D;
	    baby_step_2_r2(B, D);
	    reduce_basis(D);
	    //if(D.d0 == 5889)
	    //  cout<<"3: "<<D.d0<<" "<<D.d1<<" "<<D.d2<<endl;
	    insert(D, baby1, baby2, baby_deg, hashsize);
	  }
	  A = B;
	  reduce_basis(A);
	  //if(i==1)
	  //  cout<<"4: "<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
	} while(A.d2 < 0);
      }
    }
  }
  //GS = A;
  tempnum = A.d0;
  do {
    below(tempnum, wind, GS);
    tempnum--;
  } while(GS.d2 != wind);
}

// Baby Step-Giant Step algorithm for class number and S-regulator
// computation in unit rank 1 infrastructure.
// Input: The function field
// Output: Prints h, |h-E|, |h-E|/U, R^S, and h_x, along with 
// timing data and a factorization of h.
// See Alg. 6.2.6 of [L09].

void bsgs1(){
  ZZ R, order;
  infrastructure_ideal A, B, I, GS;
  time_t set, bstime, gstime, checktime;
  long place1, place2, k;
  int found = 0;
  ZZ_pX s, s1, s2;

  double n = tau2();
  cubic_ideal T;

  L = CeilToZZ(to_RR(L)*sqrt(n));

  reduce_basis(A);

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

  // Initialize the hashtables.
  for(long i=0; i<hashsize; i++){
    baby_deg[i] = N1;
    clear(baby1[i]);
    clear(baby2[i]);
  }

  cout<<"Number of baby steps = "<<L<<endl;

  cout<<"Making Baby Steps ... \n";
  set = time(NULL);
  order = makeBabySteps1(baby1, baby2, baby_deg, hashsize, GS);
  bstime = time(NULL);

  if(!IsZero(order)){
    if(order < (E-(genus+2))){ 
      gstime = time(NULL);
      cout<<"\n"<<"The x-Regulator is R_x = "<<2*order<<".\n";
      cout<<"The S-Regulator is R^S = "<<order<<".\n\n";
    }
    else{ 
      gstime = time(NULL);
      goto end;
    }

    //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
    cout<<" |h-E|  = "<<abs(order-E)<<endl; 
    //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
    //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
    cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<endl;
    cout<<"\n"<<"Baby Steps completed in "<<bstime-set<<" seconds.\n"<<endl;

    return;
  }
  cout<<"\n"<<"Computed "<<L<<" Baby Steps in "<<bstime-set<<" seconds. \nMaking Giant Steps ... \n";

  A = infrastructure_ideal(GS);
  k = 0;

  do{
    place1 = search(A, baby1, baby2, baby_deg, hashsize);
    if(place1 >= 0){
      found = 1;
    }
    
    // Search the hash table for the inverse of z.
    inverse(I, A);
    place2 = search(I, baby1, baby2, baby_deg, hashsize);
    if(place2 >= 0){
      found = 2;
    }
    
    B = infrastructure_ideal(A);
    A = B*GS;

    k++;
      
  }while((place1 < 0) && (place2 < 0));   
  gstime = time(NULL); 
  cout<<"Computed "<<k<<" Giant Steps in "<<gstime-bstime<<" seconds."<<endl;
  
  if(found == 1)
    order = (baby_deg[place1] - B.d0)/2;
  else{
    order = (baby_deg[place2] - I.d0)/2;
  }

 end:

  cout<<"\n"<<"A multiple of R^S is h_0 = "<<order<<".\n"; 
  cout<<"Phase 3: Completed in "<<gstime-set<<" seconds.\n"<<endl;
  
  //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<" |h-E|  = "<<abs(order-E)<<endl; 
  //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<"\n"<<endl;
  
  // Extract R^S as a factor of h_0.
  vec_ZZ factors;
  vec_long exponents;
  int i;

  R = extract(order, GS.d0/2, factors, exponents, manext);

  checktime = time(NULL);

  cout<<endl;
  cout<<"h_0     = "<<order<<" = ";
  for(i=0; i < factors.length()-1; i++){
    if(exponents[i] == 1)
      cout<<factors[i]<<" * ";
    else
      cout<<factors[i]<<"^"<<exponents[i]<<" * ";
  }
  if(exponents[i] == 1)
    cout<<factors[i];
  else
    cout<<factors[i]<<"^"<<exponents[i];
  if(!ProbPrime(factors[i]))
    cout<<" - Last factor not prime"<<endl;
  else
    cout<<endl;

  cout<<"R^S     = "<<R<<endl;
  cout<<"h^*     = "<<order/R<<endl;
  cout<<endl;

  cout<<"Phase 4: Completed in "<<checktime-gstime<<" seconds."<<endl;

  return;
}

// Baby Step-Giant Step algorithm for class number and regulator
// computation in unit rank 2 infrastructure.
// Input: The function field
// Output: Prints h, |h-E|, |h-E|/U, R, and h_x, along with 
// timing data and a factorization of h.
// A generalization of Alg. 6.2.6 of [L09].
void bsgs2(ZZ &a, ZZ &l){
  ZZ R, order;
  infrastructure_ideal A, AI, BI, B, iI, iGS, C, D;
  cubic_ideal I, GS;
  ZZ Id0, Id1, GSd0, GSd1;
  time_t set, bstime, gstime, checktime;
  long place1, place2, k;
  int found = 0;
  ZZ_pX s, s1, s2;
  ZZ tempnum, wind;
  long wind0, wind1, wind2;

  headwind2(wind0, wind1, wind2);
  wind = to_ZZ(wind2);

  double n = tau4();
  cubic_ideal T;

  L = CeilToZZ(to_RR(L)*sqrt(n));

  reduce_basis(A);
  reduce_basis(AI);

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

  // Initialize the hashtables.
  for(long i=0; i<hashsize; i++){
    baby_deg[i] = N1;
    clear(baby1[i]);
    clear(baby2[i]);
  }

  cout<<"Number of baby steps = "<<L<<endl;

  cout<<"Making Baby Steps ... \n";
  set = time(NULL);
  order = makeBabySteps2(baby1, baby2, baby_deg, hashsize, iGS);
  bstime = time(NULL);
  inf_to_cubic(iGS, GS);
  GSd0 = iGS.d0;
  GSd1 = iGS.d1;

  if(!IsZero(order)){
    if(order < (E-genus)){ 
      gstime = time(NULL);
    }
    else{ 
      gstime = time(NULL);
      goto end;
    }

    //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
    cout<<" |h-E|  = "<<abs(order-E)<<endl; 
    //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
    //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
    cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<endl;
    cout<<"\n"<<"Baby Steps completed in "<<bstime-set<<" seconds.\n"<<endl;

    return;
  }
  cout<<"\n"<<"Computed "<<L<<" Baby Steps in "<<bstime-set<<" seconds. \nMaking Giant Steps ... \n";

 preloop:

  A = infrastructure_ideal(iGS);
  k = 0;
  if(A.d2 > 0){
    reduce_basis1(A);
    while(A.d2 >= 0){
      C=A;
      baby_step_1_r2(C, A);
    }
    reduce_basis(C);
    A = C;
  }
  AI = A;
  // Set the inverse giant step.
  C = A;
  tempnum = C.d0;
  do{
    inverse(iI, C);
    if(iI.d2 > wind){
      reduce_basis1(iI);
      while(iI.d2 >= wind){
	C=iI;
	baby_step_1_r2(C, iI);
      }
      reduce_basis(C);
      iI = C;
    }
    else if(iI.d2 < wind){
      reduce_basis2(iI);
      while(iI.d2 <= wind){
	C=iI;	
	baby_step_2_r2(C, iI);
      }
      reduce_basis(C);
      iI = C;
      if(iI.d2 < wind){
	do{
	  baby_step_0_r2(iI, C);
	  iI = C;
	  reduce_basis2(C);
	  while(C.d2 <= wind){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  iI = D;
	  reduce_basis(iI);
	} while(iI.d2 < wind);
      }
    }
    else;
    tempnum--;
    below(tempnum, Zero, C);
  } while(iI.d2 != wind);

  inf_to_cubic(iI, I);
  Id0 = iI.d0;
  Id1 = iI.d1;

 loop:
  do{
    //if(k == 187){
    //  cout<<k<<"-> "<<A.d0<<" "<<A.d1<<" "<<A.d2<<" "<<A.d<<endl;
    //  A.print();
    //}
    place1 = search(A, baby1, baby2, baby_deg, hashsize);
    if(place1 >= 0){
    check1:
      order = baby_deg[place1] - A.d0;
      //cout<<order<<" "<<check_inf2(order)<<" "<<checkPrimes(order)<<" "<<order%l<<endl;
      if(order%l == a){
	if(check_inf2(order)){
	  if(checkPrimes(order))
	    found = 1;
	  else{
	    place1 = -1;
	    order = 0;
	  }
	}
	else{
	  place1 = -1;
	  order = 0;
	}
      }
      else{
	place1 = -1;
	order = 0;
      }
    }
    
    BI = infrastructure_ideal(AI);
    
    giant_step_r2(AI, BI, I, Id0, Id1, wind);
    //AI = BI*I;
    if(AI.d2 < 0){     
      while(AI.d2 <= 0){
	reduce_basis2(AI);
	C=AI;	
	baby_step_2_r2(C, AI);
	reduce_basis(AI);
	place2 = search(AI, baby1, baby2, baby_deg, hashsize);
	if(place2 >= 0)
	  goto check2;
      }
      reduce_basis(C);
      AI = C;
      if(AI.d2 < 0){
	do{
	  baby_step_0_r2(AI, C);
	  AI = C;
	  place2 = search(AI, baby1, baby2, baby_deg, hashsize);
	  if(place2 >= 0)
	    goto check2;
	  reduce_basis2(C);
	  while(C.d2 <= 0){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  AI = D;
	  reduce_basis(AI);
	  place2 = search(AI, baby1, baby2, baby_deg, hashsize);
	  if(place2 >= 0)
	    goto check2;
	} while(AI.d2 < 0);
      }
    }
    else if(AI.d2 > 0){
      reduce_basis1(AI);
      while(AI.d2 >= 0){
	C=AI;
	baby_step_1_r2(C, AI);
      }
      reduce_basis(C);
      AI = C;
    }
    //cout<<k<<"-> "<<AI.d0<<" "<<AI.d1<<" "<<AI.d2<<" "<<A.d<<endl;
    place2 = search(AI, baby1, baby2, baby_deg, hashsize);
    if(place2 >= 0){
    check2:
      order = baby_deg[place2] - AI.d0;
      //cout<<order<<" "<<check_inf2(order)<<" "<<checkPrimes(order)<<" "<<order%l<<endl;
      if(order%l == a){
	if(check_inf2(order)){
	  if(checkPrimes(order))
	    found = 2;
	  else{
	    place2 = -1;
	    order = 0;
	  }
	}
	else{
	  place2 = -1;
	  order = 0;
	}
      }
      else{
	place2 = -1;
	order = 0;
      }
    }
    
    B = infrastructure_ideal(A);
    giant_step_r2(A, B, GS, GSd0, GSd1, wind);
    //A = B*GS;

    if(A.d2 < 0){
      while(A.d2 <= 0){
	reduce_basis2(A);
	C=A;	
	baby_step_2_r2(C, A);
	reduce_basis(A);
	place1 = search(A, baby1, baby2, baby_deg, hashsize);
	if(place1 >= 0){
	  k++;
	  goto check1;
	}
      }
      reduce_basis(C);
      A = C;
      if(A.d2 < 0){
	do{
	  baby_step_0_r2(A, C);
	  A = C;
	  place1 = search(A, baby1, baby2, baby_deg, hashsize);
	  if(place1 >= 0){
	    k++;
	    goto check1;
	  }
	  reduce_basis2(C);
	  while(C.d2 <= 0){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  A = D;
	  reduce_basis(A);
	  place1 = search(A, baby1, baby2, baby_deg, hashsize);
	  if(place1 >= 0){
	    k++;
	    goto check1;
	  }
	} while(A.d2 < 0);
      }
    }
    else if(A.d2 > 0){
      reduce_basis1(A);
      while(A.d2 >= 0){
	C=A;
	baby_step_1_r2(C, A);
      }
      reduce_basis(C);
      A = C;
    }
    k++;
      
  }while((place1 < 0) && (place2 < 0) && (A.d0 < U)); 
  if(A.d0 >= U){      
    cout<<"Making new giant step."<<endl;
    GSd0 -= genus;
    do { 
      below(GSd0, wind, iGS);
      GSd0--;
    } while(iGS.d2 != wind);
    inf_to_cubic(iGS, GS);
    GSd0 = iGS.d0;
    GSd1 = iGS.d1;
    goto preloop;
  }
  

  gstime = time(NULL); 
  cout<<"Computed "<<k<<" Giant Steps in "<<gstime-bstime<<" seconds."<<endl;

  //if(found == 1)
  //  order = baby_deg[place1] - B.d0;
  //else{
  //  order = baby_deg[place2] - AI.d0;
  //}

 end:

  //cout<<"Order OK? "<<check_inf2(order)<<endl;

  cout<<"\n"<<"A multiple of R is h = "<<order<<".\n"; 
  cout<<"Phase 3: Completed in "<<gstime-set<<" seconds.\n"<<endl;
  
  //cout<<"   |h-E1|   = "<<abs(order-E1)<<endl;
  cout<<" |h-E|  = "<<abs(order-E)<<endl; 
  //cout<<"|h-E1|/L1^2 = "<<to_RR(abs(order-E1))/to_RR(L12)<<endl;
  //cout<<"|h-E2|/L2^2 = "<<to_RR(abs(order-E))/to_RR(L2)<<endl;
  cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<"\n"<<endl;

  // Extract a multiple of R as a factor of h_0.
  vec_ZZ factors;
  vec_long exponents;
  int i;

  R = extract2(order, GSd0, factors, exponents, manext);

  // If R == order, then we know that h_x = 1 and we are done.
  // Otherwise, the regulator is some multiple of R. Determine which here.
  // Need to add in functionality to run the Kangaroo method if 
  // there is insufficient memory.
  if(R != order){
    ZZ icn = order/R;
    cout<<"1 <= h_x <= "<<icn<<". Searching the infrastructure to determine h_x and R."<<endl;
    if(!latticeSearch(icn, R)){
      place1 = place2 = -1;
      goto loop;
    }
  }
  else{
    ZZ icn = order/R;
    latticeSearch(icn, R);
  }
  checktime = time(NULL);

  cout<<"h     = "<<order<<" = ";
  if(factors.length() == 1){
    cout<<"prime"<<endl;
  }
  else{
    for(i=0; i < factors.length()-1; i++){
      if(exponents[i] == 1)
	cout<<factors[i]<<" * ";
      else
	cout<<factors[i]<<"^"<<exponents[i]<<" * ";
    }
    if(exponents[i] == 1)
      cout<<factors[i];
    else
      cout<<factors[i]<<"^"<<exponents[i];
    if(!ProbPrime(factors[i]))
      cout<<" - Last factor not prime"<<endl;
    else
      cout<<endl;
  }
  cout<<"R     = "<<R<<endl;
  cout<<"h_x   = "<<order/R<<endl;
  cout<<endl;

  cout<<"Phase 4: Completed in "<<checktime-gstime<<" seconds."<<endl;

  return;
}

// latticeSearch() applies the BSGS method in the Jacobian lattice
// with 2-distance of 1, applying more baby steps up to a 2-distance
// of icn. When a unit (the identity) is found, we can quickly find
// the actual ideal class number and regulator.
// Input: icn - potential ideal class number
//            - the 2-distance of a known unit
//        R - potential regulator
//          - the exponent of the Jacobian
// Output: The actual ideal class number, icn, and regulator, R, 
//         or an error if there's not enough memory.
// This function needs fixing for the case in which there's not enough 
// memory to perform baby steps.
// May need to modify it to optimize via [F09].

int latticeSearch(ZZ &icn, ZZ &R){
  ZZ bsNum = RoundToZZ(sqrt((tau5()+to_RR(icn))*to_RR(R)/2.0));
  ZZ d2factor;
  ZZ unit2;
  long hashsize = NextPrime(4*to_long(bsNum));
  infrastructure_ideal A, B, C, iGS, D;
  cubic_ideal GS;
  ZZ GSd0, GSd1, wind;
  int k=0, place1=-1, place2=-1, found = 0;
  int repeat = 1;
  //hashbits = NumBits(hashsize);

  // If there is enough room to use BSGS, do so.
  // We will only use at most a third of available memory.
  if((4*bsNum*(2*q.size() + R.size())*NTL_ZZ_NBITS/24 < 1048576*memory) && (icn < 10000)){
      
    ZZVec baby1 = ZZVec(hashsize, q.size());  // First and
    ZZVec baby2 = ZZVec(hashsize, q.size());  // second baby step arrays
    ZZVec baby_deg = ZZVec(hashsize, R.size()); // Stores the degrees.
    if((&baby1 == NULL) || (&baby2 == NULL) || (&baby_deg == NULL) ){
      cout<<"Not enough memory available to create baby steps.\n";
      cout<<"Hash table size = "<<hashsize<<".\n"<<"Total memory required = "<<3*hashsize*sizeof(ZZ)<<" bytes.\n\n";
      cout<<"Using Kangaroo Method.\n\n";
      latticeRoo(icn, R);
    }
    
    // Initialize the hashtables.
    for(long i=0; i<hashsize; i++){
      baby_deg[i] = N1;
      clear(baby1[i]);
      clear(baby2[i]);
    }

    cout<<"Making "<<bsNum<<" Baby Steps ... \n";
    //set = time(NULL);
    latticeBS(bsNum, baby1, baby2, baby_deg, hashsize, iGS);
    //bstime = time(NULL);

    inf_to_cubic(iGS, GS);
    GSd0 = iGS.d0;
    GSd1 = iGS.d1;
    wind = iGS.d2;
    cout<<"\nCompleted baby steps. Making giant steps."<<endl;

  start:
    
    do{
      reduce_basis(A);
      baby_step_0_r2(A, B);
      reduce_basis2(B);
      baby_step_2_r2(B, A);
    } while (A.d2 < 1);

    k = 0;
    place1 = search(A, baby1, baby2, baby_deg, hashsize);
    if(place1 >= 0){
      found = 1;
    }
    
    while((place1 < 0) && (place2 < 0) && (A.d0 <= R)){
      k++;
      
      B = infrastructure_ideal(A);
      giant_step_r2(A, B, GS, GSd0, GSd1, wind);
      //A = B*GS;
          
      if(A.d2 < 1){
	while(A.d2 <= 1){
	  reduce_basis2(A);
	  C=A;	
	  baby_step_2_r2(C, A);
	  reduce_basis(A);
	  place1 = search(A, baby1, baby2, baby_deg, hashsize);
	  if(place1 >= 0)
	    goto check1;
	}
	reduce_basis(C);
	A = C;
	if(A.d2 < 1){
	  do{
	    baby_step_0_r2(A, C);
	    A = C;
	    place1 = search(A, baby1, baby2, baby_deg, hashsize);
	    if(place1 >= 0)
	      goto check1;
	    reduce_basis2(C);
	    while(C.d2 <= 1){
	      D = C;
	      baby_step_2_r2(D, C);
	    }
	    A = D;
	    reduce_basis(A);
	    place1 = search(A, baby1, baby2, baby_deg, hashsize);
	    if(place1 >= 0)
	      goto check1;
	  } while(A.d2 < 1);
	}
      }
      else if(A.d2 > 1){
	reduce_basis1(A);
	while(A.d2 >= 1){
	  C=A;
	  baby_step_1_r2(C, A);
	}
	reduce_basis(C);
	A = C;
      }
      //if((k>900) && (k < 930))
      //	cout<<k<<": "<<A.d0<<" "<<A.d0+GSd0<<" "<<A.d2<<" "<<A.d<<endl;
      place1 = search(A, baby1, baby2, baby_deg, hashsize);
      
      // Search vertically until h_x.
      B = A;
      while((B.d2 < icn) && (place2 < 0)){
	reduce_basis2(B);
	baby_step_2_r2(B, C);
	B = C;
	reduce_basis(C);
	place2 = search(C, baby1, baby2, baby_deg, hashsize);
      }
    }
   
    cout<<endl;
    
    //cout<<A.d0<<" "<<R<<" "<<place1<<" "<<place2<<endl;
    //if(place2 != -1)
    // cout<<C.d0<<" "<<baby_deg[place2]<<" "<<C.d2<<endl;
    if(place1 >= 0){
    check1:
      unit2 = A.d0 - baby_deg[place1];
      orthogonalize(R, unit2, One);
      printUnits(R, Zero, unit2, One);
      return 1;
    }
    else if(place2 >= 0){
      unit2 = C.d0 - baby_deg[place2];
      // Just to make sure that we get the 2-distance of the unit correct.
      // There's an extremely small chance that A.d2 != C.d2.
      //cout<<unit2<<" "<<C.d2<<endl;
      below(unit2, C.d2, A);
      //cout<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
      if(IsOne(A.d))
	cout<<"Unit found at "<<A.d0<<"(oo_1 - oo_0) + "<<A.d2<<"(oo_1 - oo_2).\n\n";
      else{
	int l = 1;
	while(!IsOne(A.d)){
	  C.d2 += l;
	  below(unit2, C.d2, A);
	  //cout<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
	  if(l>0) { l++;  l = -l; }
	  else { l--; l = -l; }
	}
	cout<<"Unit found at "<<A.d0<<"(oo_1 - oo_0) + "<<A.d2<<"(oo_1 - oo_2).\n\n";
      }
      
      // Given the unit that we've found, we search the lattice
      // for the unit with minimal positive 2-distance in order
      // to nail down the regulator and ideal class number exactly.
      if(processUnit(icn, R, unit2, C.d2))
	return 1;
      else
	return 0;
    }
    else{
      if(repeat){
	cout<<"Restarting Giant Steps with new giant step."<<endl;
	ZZ tempnum = iGS.d0 - 2*genus;
	do {
	  below(tempnum, wind, iGS);
	  tempnum--;
	} while(iGS.d2 != wind);
	inf_to_cubic(iGS, GS);
	GSd0 = iGS.d0;
	GSd1 = iGS.d1;
	wind = iGS.d2;
	A = infrastructure_ideal();
	repeat = 0;
	goto start;
      }
      else{
	if(!latticeRoo(icn, R))
	  return 0;
      }
      //R *= icn;   icn = 1;    cout<<endl;
    }
  }
  else{
    if(!latticeRoo(icn, R))
      return 0;
  }

  return 1;
}

// latticeRoo() applies the Kangaroo method in the Jacobian lattice
// with 2-distance of icn. When a unit (the identity) is found, we can 
// quickly find the actual ideal class number and regulator.
// Input: icn - potential ideal class number
//            - the 2-distance of a known unit
//        R - potential regulator
//          - the exponent of the Jacobian
// Output: The actual ideal class number, icn, and regulator, R, 
// This function is not complete.

int latticeRoo(ZZ &icn, ZZ &R){
  int i, found = 0;
  int distbits; // The number of zero bits to determine 
                // if an ideal is distinguished.
  int distbits2;
  int dbadd = 4;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  ZZ avgjump;  // The average of the jumps.
  ZZ randjump; // The upper bound on the random jump distance.
  infrastructure_ideal tame, wild;
  infrastructure_ideal temp, C, D;
  cubic_ideal jumps[64];
  ZZ jumpdistance[64];
  ZZ jumpdistance1[64];
  ZZ jumpsum = ZZ();

  ZZ matchdist = ZZ();

  ZZ bsjumptotal = ZZ(); // The number of baby step  jumps for each kangaroo.
  ZZ gsjumptotal = ZZ(); // The number of giant step jumps for each kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  long expdps; // The expected number of distinguished points.
  long hashsize;
  ZZ order; // The order of the group.
  double tauratio = tau5(); // T_G/T_B - giant step to baby step ratio.
  time_t set, rootime, checktime;
  ZZ Stau;
  long dist0, dist1, dist2;
  ZZ wind0, wind2; 

  headwind2(dist0, dist1, dist2);
  wind0 = to_ZZ(dist0);
  wind2 = to_ZZ(dist2);

  cout<<"Initializing ...";
  cout<<flush;

  // Initialize the kangaroos.
  reduce_basis(tame);
  reduce_basis(wild);
  baby_step_0_r2(tame, temp);
  reduce_basis2(temp);
  while(tame.d2 <= icn){
    temp = tame;
    baby_step_2_r2(temp, tame);
  }
  tame = temp;

  Stau = RoundToZZ(to_RR(q)/tauratio);

  if(q < 1200)
    dbadd-=2;
  if(q < 120)
    dbadd-=2;

  // Set the average jump distance and the jump distances, 
  avgjump = RoundToZZ(sqrt(to_RR(R/2)*(2.0*tauratio-1.0)) - tauratio + 1.0);

  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 1 + wind0);
    
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*(avgjump)) || (jumpsum > (i+1)*(avgjump)) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 1 + wind0);
      jumpsum += jumpdistance[i];
    }
    
    jumpdistance[63] -= jumpdistance[i];
  }
  jumpdistance[63] += (64*(avgjump) + wind0 + 32);

  // Set the number of bits for a distinguished ideal.
  distbits = (NumBits(avgjump)-1)/2 - 2 + dbadd;

  distbits += (distbits%2);
  distbits2 = distbits/2;

  RR expgsjumps = 4.0*sqrt(to_RR(R/2)/(2.0*tauratio-1.0)) + to_RR(2*power2_ZZ((long)distbits))/tauratio;

  // Set the hashtables and their size.
  hashsize = NextPrime(to_long(10*RightShift(RoundToZZ(expgsjumps), distbits)));
  hashbits = NumBits(hashsize);

  // The hash tables for the tame kangaroo.

  ZZVec dpt1 = ZZVec(hashsize, q.size());  // store s(0)
  ZZVec dpt2 = ZZVec(hashsize, q.size());  // store s(1)
  ZZVec dptd = ZZVec(hashsize, R.size());  // store the distance

  // The hash tables for the wild kangaroo.

  ZZVec dpw1 = ZZVec(hashsize, q.size());  // store s(0)
  ZZVec dpw2 = ZZVec(hashsize, q.size());  // store s(1)
  ZZVec dpwd = ZZVec(hashsize, R.size());  // store the distance

  // Initialize the hashtables.
  for(long j=0; j<hashsize; j++){
    dptd[j] = dpwd[j] = N1;
  }
  // Initialize the jumps.
  for(i=0; i<64; i++){
    below(jumpdistance[i], wind2, temp);
    while(temp.d2 != wind2){
      jumpdistance[i]--;
      below(jumpdistance[i], wind2, temp);
    }
    jumpdistance1[i] = temp.d1;
    inf_to_cubic(temp, jumps[i]);
  }

  cout<<" Done. Beginning the jumping."<<endl;
  //cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  // cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  //cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  expdps = to_long(RightShift(RoundToZZ(expgsjumps), distbits));
  cout<<"Distinguished points expected:   "<<expdps<<endl;
  
  set = time(NULL);
  //cout<<"Initializations completed in "<<set-timingtest<<" seconds"<<endl;
  while(!found){
    // FOR THE TAME KANGAROO:
    /************************************************/
    /* Make baby steps until L(tame)(0) < [q/tau_3]. */
    /************************************************/
    while(rep(coeff(tame.d, 0)) >= Stau){
      baby_step_0_r2(tame, temp);
      tame = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to icn as possible.
      if(tame.d2 < icn){
	reduce_basis2(tame);
	while(tame.d2 <= icn){
	  C=tame;	
	  baby_step_2_r2(C, tame);
	}
	reduce_basis(C);
	tame = C;
	if(tame.d2 < icn){
	  do{
	    baby_step_0_r2(tame, C);
	    tame = C;
	    reduce_basis2(C);
	    while(C.d2 <= icn){
	      D = C;
	      baby_step_2_r2(D, C);
	    }
	    tame = D;
	    reduce_basis(tame);
	  } while(tame.d2 < icn);
	}
      } 
      bsjumptotal++;
    }
    
    // I'll only check for DPs in the set S_tau.
    // Check if the ideal is distinguished.
    roocoeff = rep(coeff(tame.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(tame.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(tame, 0, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, matchdist, hashsize); 
	distpoints++;
        
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(tame.d0-matchdist);
	  below(order, icn, C);
	  if(IsOne(C.d)){
	    cout<<endl<<endl;
	    cout<<"Unit found at "<<order<<"(oo_1 - oo_0) + "<<icn<<"(oo_1 - oo_2).\n\n";
	    processUnit(icn, R, order, icn);
	    return 1;
	  }
	  else
	    found = 0; // It's a dud.
	}
	if(distpoints > 10*expdps)
	  return 0;
      }
    }
        
    // FOR THE WILD KANGAROO:
    // Make baby steps until s(0) = 0 (mod tau).
    /************************************************/
    /* Make baby steps until L(wild)(0) < [q/tau_3]. */
    /************************************************/
    while(rep(coeff(wild.d, 0)) >= Stau){
      baby_step_0_r2(wild, temp);
      wild = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to icn as possible.
      if(wild.d2 < 0){
	reduce_basis2(wild);
	while(wild.d2 <= 0){
	  C=wild;	
	  baby_step_2_r2(C, wild);
	}
	reduce_basis(C);
	wild = C;
	if(wild.d2 < 0){
	  do{
	    baby_step_0_r2(wild, C);
	    wild = C;
	    reduce_basis2(C);
	    while(C.d2 <= 0){
	      D = C;
	      baby_step_2_r2(D, C);
	    }
	    wild = D;
	    reduce_basis(wild);
	  } while(wild.d2 < 0);
	}
      }
      bsjumptotal++;
    }
    
    // Check if the ideal is distinguished.
    roocoeff = rep(coeff(wild.mu0, 0));
    trunc(dptest, roocoeff, distbits2);
    if(IsZero(dptest)){
      roocoeff = rep(coeff(wild.nu0, 0));
      trunc(dptest, roocoeff, distbits2);
      if(IsZero(dptest)){
	found = distpoint(wild, 1, dpt1, dpt2, dptd, dpw1, dpw2, dpwd, matchdist, hashsize); 
	distpoints++;
	
	printf("\r");
	printf("Distinguished points: %ld", distpoints);
	cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(wild.d0-matchdist);
	  below(order, icn, C);
	  if(IsOne(C.d)){
	    cout<<endl<<endl;
	    cout<<"Unit found at "<<order<<"(oo_1 - oo_0) + "<<icn<<"(oo_1 - oo_2).\n\n";
	    processUnit(icn, R, order, icn);
	    return 1;
	  }
	  else
	    found = 0; // It's a dud.
	}
	if(distpoints > 10*expdps)
	  return 0;
      }
    }
    
    // Make the jumps.
    i = vmap(tame);
    temp = tame;
    giant_step_r2(tame, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    //tame = temp*jumps[i];

    i = vmap(wild);
    temp = wild;
    giant_step_r2(wild, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    //wild = temp*jumps[i];

    // Correct to make sure that the 2-component of distance is 
    // as close to icn as possible.
    if(tame.d2 < icn){
      reduce_basis2(tame);
      while(tame.d2 <= icn){
	C=tame;	
	baby_step_2_r2(C, tame);
      }
      reduce_basis(C);
      tame = C;
      if(tame.d2 < icn){
	do{
	  baby_step_0_r2(tame, C);
	  tame = C;
	  reduce_basis2(C);
	  while(C.d2 <= icn){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  tame = D;
	  reduce_basis(tame);
	} while(tame.d2 < icn);
      }
    }
    else if(tame.d2 > icn){
      reduce_basis1(tame);
      while(tame.d2 >= icn){
	D=tame;
	baby_step_1_r2(D, tame);
      }
      reduce_basis(D);
      tame = D;
    }

    if(wild.d2 < 0){
      reduce_basis2(wild);
      while(wild.d2 <= 0){
	C=wild;	
	baby_step_2_r2(C, wild);
      }
      reduce_basis(C);
      wild = C;
      if(wild.d2 < 0){
	do{
	  baby_step_0_r2(wild, C);
	  wild = C;
	  reduce_basis2(C);
	  while(C.d2 <= 0){
	    D = C;
	    baby_step_2_r2(D, C);
	  }
	  wild = D;
	  reduce_basis(wild);
	} while(wild.d2 < 0);
      }
    }
    else if(wild.d2 > icn){
      reduce_basis1(wild);
      while(wild.d2 >= icn){
	D=wild;
	baby_step_1_r2(D, wild);
      }
      reduce_basis(D);
      wild = D;
    }

    gsjumptotal+=2;     
  }
}
// Given a unit with a nonzero 2-distance, processUnit() searches the 
// Jacobian lattice for the second unit in a fundamental system.
// This will identify the ideal class number and regulator precisely.
// Input: icn - potential ideal class number; h/exponent
//        R - potential regulator; the exponent of hte Jacobian.
//        unit2 - the 0-distance of a found unit.
//        dist2 - the 2-distance of a found unit.

int processUnit(ZZ &icn, ZZ &R, ZZ &unit2, ZZ &dist2){
  int i, j; 
  ZZ d2factor;
  infrastructure_ideal C;
  for(j=1; j<dist2; j++){
    if(((icn%j) == 0) && ((dist2%j) == 0)){
      for(i = 0; i < dist2/j; i++){
	if(IsZero((unit2+i*R)%(dist2/j))){
	  unit2 += i*R; unit2 /= (dist2/j); 
	  //cout<<i<<" "<<unit2<<" "<<j<<endl;
	  d2factor = to_ZZ(j); 
	  below(unit2, d2factor, C);
	  //cout<<C.d0<<" "<<C.d1<<" "<<C.d2<<" "<<C.d<<endl;
	  if(IsOne(C.d)){
	    orthogonalize(R, C.d0, C.d2);
	    printUnits(R, Zero, C.d0, C.d2);
	    R *= C.d2;   icn /= C.d2; 
	    return 1;
	  }
	  unit2 *= (dist2/j); unit2 -= i*R;
	}
      }
    }
  }
  if(!IsZero(icn%dist2))
    return 0;

  orthogonalize(R, unit2, dist2);
  printUnits(R, Zero, unit2, dist2);
  R*=dist2;  icn /= dist2;  
  
  return 1;
}

// Given a system of fundamental units, finds a system as orthogonal as possible.

void orthogonalize(ZZ &e10, ZZ &e20, ZZ &e22){
  mat_ZZ unitMatrix;
  ZZ det2 = sqr(e10*e22);

  unitMatrix.SetDims(2,2);
  unitMatrix[0][0] = e10;
  unitMatrix[0][1] = Zero;
  unitMatrix[1][0] = e20;
  unitMatrix[1][1] = e22;

  LLL(det2, unitMatrix);

  u10 = unitMatrix[0][0];
  u12 = unitMatrix[0][1];
  u20 = unitMatrix[1][0];
  u22 = unitMatrix[1][1];
}

// Prints the units of the function field via their divisors.
// Input: e10 - 0-distance of the first unit.
//        e12 - 2-distance of the first unit.
//        e20 - 0-distance of the second unit.
//        e22 - 2-distance of the second unit.
void printUnits(ZZ &e10, ZZ &e12, ZZ &e20, ZZ &e22){
  // Print the first unit.
  if(IsZero(e12))
    cout<<"Units: div(unit_1) = "<<e10<<"(oo_1 - oo_0)"<<endl;
  else
    cout<<"Units: div(unit_1) = "<<e10<<"(oo_1 - oo_0) + "<<e12<<"(oo_1 - oo_2)"<<endl;
  // Print the second unit.
  if(IsOne(e22))
    cout<<"       div(unit_2) = "<<e20<<"(oo_1 - oo_0) + (oo_1 - oo_2)"<<endl;
  else
    cout<<"       div(unit_2) = "<<e20<<"(oo_1 - oo_0) + "<<e22<<"(oo_1 - oo_2)"<<endl;
  cout<<endl;
  cout<<"Orthogonalized:"<<endl;
  cout<<"Units: div(unit_1) = "<<u10<<"(oo_1 - oo_0) + "<<u12<<"(oo_1 - oo_2)"<<endl;
  cout<<"       div(unit_2) = "<<u20<<"(oo_1 - oo_0) + "<<u22<<"(oo_1 - oo_2)"<<endl;
  cout<<endl;
}

// getinput() reads in the data from the input file, initializes a number 
// of variables for the desired computation, and constructs a polynomial 
// to generate the function field if needed.
// Input: input - the file that we are reading data from.
int getinput(char *input){
  FILE *inputptr;
  char thisLine[128];
  char token[64], value[64];
  int dG=4, dH=0, rand=2;
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
    else if(strncmp(token, "memory:", 7)==0)
      memory = atoi(value);
    else if(strncmp(token, "extract:", 8)==0)
      manext = atoi(value);
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
      G = ZZ_pX();
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
      cout<<"Generating new G(x). Original was either not squarefree or not coprime to H(x).\n\n";
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
  for(i=0; i<21; i++)
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
    else if(isCube(LeadCoeff(f)))
      rank = 1;
    else
      rank = 0;
  }

  alpha = alphahat();
  
  fclose(inputptr);
 
  // For hash storage. 

  if(genus < 5)
    type = 1;
  else
    type = 2;

  return(1);
 
}

// printInitData() prints information about the function field that
//  we are working with: the polynomial, ground field, unit rank, and
//  the genus.

void printInitData(){
  int dD = deg(f), dG = deg(G), dH = deg(H);
  int i;

  cout<<"\n";
  cout<<"Working with the function field F_"<<q<<"(x, y), where"<<endl;
  cout<<"y^3 = G(x) H^2(x) = ";

  if (dD > 1){
    if(!IsOne(coeff(f,dD)))
      cout<<coeff(f,dD);
    cout<<"x^"<<dD;
  }
  for(i=dD-1; i>1; i--){
    if(!IsZero(coeff(f, i))) {
      cout<<" + ";
      if(!IsOne(coeff(f,i)))
	cout<<coeff(f,i);
      cout<<"x^"<<i;
    }
  }
  if(!IsZero(coeff(f,1))){
    if(dD > 1) { cout<<" + "; }
    if(!IsOne(coeff(f,1)))
      cout<<coeff(f,1);
    cout<<"x";
  }
  if(!IsZero(coeff(f,0)))
    cout<<" + "<<coeff(f,0)<<endl;
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

// Needs to be fixed so that this file indeed performs the same task as
//  printInitData(), but printing it to a file. Currently, it prints to
//  the screen.

void printInitDataFile(char *file){
  //FILE *outputptr;
  int dD = deg(f), dG = deg(G), dH = deg(H);
  int i;

  // For now, all this prints out to the screen.
  cout<<"Working with the function field F_"<<q<<"(x, y), where"<<endl;
  cout<<"y^3 = G(x) H^2(x) = ";

  if (dD > 1){
    if(!IsOne(coeff(f,dD)))
      cout<<coeff(f,dD);
    cout<<"x^"<<dD;
  }
  for(i=dD-1; i>1; i--){
    if(!IsZero(coeff(f, i))) {
      cout<<" + ";
      if(!IsOne(coeff(f,i)))
	cout<<coeff(f,i);
      cout<<"x^"<<i;
    }
  }
  if(!IsZero(coeff(f,1))){
    if(dD > 1) { cout<<" + "; }
    if(!IsOne(coeff(f,1)))
      cout<<coeff(f,1);
    cout<<"x";
  }
  if(!IsZero(coeff(f,0)))
    cout<<" + "<<coeff(f,0)<<endl;
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

// main() has the input file read in and chooses which algorithm with
//  which to run the class number and/or regulator compuation for the 
//  given purely cubic function field.

int main(int argc, char *argv[]){
  ZZ l, a;
  time_t set, time1, time2, time3;

  if(!((argc == 2) || (argc == 3))){
    cout<<"Usage: ./classnumber inputfile [outputfile]\n\n";
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
  
  // Step 1: Compute an approximation E of the class number h
  //         and an integer L such that |h - E| < L^2.
  set = time(NULL);

  approxh(1);

  time1 = time(NULL);

  cout<<"Phase 1: Completed in "<<time1-set<<" seconds.\n\n";

  // Step 2: Use extra information about h in (E-L^2, E+L^2),
  //         such as its distribution in the interval
  //         or h mod r for small primes r.

  // Find elements of small order.
  //if(rank == 0){
  smallOrder0(a, l);
  cout<<"Phase 2: h = "<<a<<" (mod "<<l<<").\n\n";
  //}

  time2 = time(NULL);

  // Step 3: Find h in the interval (E-L^2, E+L^2) via
  //         A) Baby Step, Giant Step or
  //         B) Pollard's Kangaroo
  if(rank == 0){
    if(bsgs)
      bsgs0(a, l);
    if(kang)
      kangaroo0(a, l); 
    time3 = time(NULL);
    cout<<"Phase 3: Completed in "<<time3-time2<<" seconds.\n\n";
  
  }
  else if(rank == 1){
    // Get the basis, {1, rho, omega}.
    init(2*genus);
    
    if(bsgs)
      bsgs1();
    if(kang)
      kangaroo1();
  }
  else{
    // Get the basis, {1, rho, omega}.
    init(2*genus);
    if(bsgs)
      bsgs2(a, l);
    if(kang){
      if(IsOne(l))
	kangaroo2();
      else
	kangaroo2(a, l);
    }
  }
  time3 = time(NULL);
  cout<<"Total running time: "<<time3-set<<" seconds.\n\n";
}
