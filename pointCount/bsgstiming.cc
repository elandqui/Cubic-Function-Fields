/* 
   bsgstiming.cc

   This file contains functions to compute timing data for the
   ratio between a baby step and giant step (possibly plus an inverse).
   Timings will involve 10^6 steps each.
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
   For unit rank 2 - Need to incorporate Kangaroo methods for extracting
                     R from h for larger examples. (Phase 4)
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   4. Incorporate Kangaroo methods into Phases 3 and 4 - Check Fontein's
      paper for improving Phase 4.
   5. Organize this file better and possibly even split it up into several
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
  /*
  tau = (int)n;
  double tauprime = (double)tau; 
  avgjump = CeilToZZ(2.0*sqrt(alpha*to_RR(U)*(2.0*tauprime-1.0)));
  distbits = (NumBits(avgjump)-1)/2 + dbadd;
  distbits += (distbits%2);
  // Set tau depending on whichever is faster, rounding tau' up or down.
  // Calculate the predicted running time for each kangaroo for rounding up and down.
  RR timepred1 = (2.0*sqrt(alpha*to_RR(U))/sqrt(to_RR(2.0*tau-1.0)) + to_RR(power2_ZZ(distbits))/tauprime )*(1.0+(tauprime - 1.0)/n);

  tau++;
  tauprime = (double)tau;  
  avgjump = CeilToZZ(2.0*sqrt(alpha*to_RR(U)*(2.0*tauprime-1.0)));
  distbits = (NumBits(avgjump)-1)/2 + dbadd;
  distbits += (distbits%2);
  RR timepred2 = (2.0*sqrt(alpha*to_RR(U))/sqrt(to_RR(2.0*tau-1.0)) + to_RR(power2_ZZ(distbits))/tauprime )*(1.0+(tauprime - 1.0)/n);
  
  tau++;
  tauprime = (double)tau;  
  avgjump = CeilToZZ(2.0*sqrt(alpha*to_RR(U)*(2.0*tauprime-1.0)));
  distbits = (NumBits(avgjump)-1)/2 + dbadd;
  distbits += (distbits%2);
  RR timepred3 = (2.0*sqrt(alpha*to_RR(U))/sqrt(to_RR(2.0*tau-1.0)) + to_RR(power2_ZZ(distbits))/tauprime )*(1.0+(tauprime - 1.0)/n);
  

  // Adjust tau depending on what would be faster.
  if((timepred2 < timepred1) && (timepred2 < timepred3))
    tau--;
  else if(timepred1 < timepred3)
    tau-=2;
  else;
  */

  // Set the average jump distance, the jump distances, 
  // and the starting position for tame.
  avgjump = RoundToZZ(2.0*sqrt(alpha*to_RR(U)*(2.0*tauratio-1.0)) - 2.0*(tauratio-1.0));

  if(genus%3 == 1)
    headwind = (genus+2)/3;
  else
    headwind = genus/3;

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
  
  ZZ one = to_ZZ(1);

  ZZ R = extract(order, one, factors, exponents, manext);
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

  //for(i=0; i<=genus; i++)
  //  cout<<"Frequency of headwind strength "<<i<<": "<<distdiffs[i]<<endl;
  //cout<<endl;
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
  infrastructure_ideal temp, D;
  infrastructure_ideal jumps[64];
  ZZ jumpdistance[64];
  ZZ jumpsum = ZZ();
  int j, count = 0;

  long headwind;

  ZZ matchdist = ZZ();

  ZZ bsjumptotal = ZZ(); // The number of baby step  jumps for each kangaroo.
  ZZ gsjumptotal = ZZ(); // The number of giant step jumps for each kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  long hashsize;
  ZZ order; // The order of the group.
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
  /*
  tau = (int)n;
  double tauprime = (double)tau; 
  avgjump = CeilToZZ(2.0*sqrt(alpha*to_RR(U)*(2.0*tauprime-1.0)));
  distbits = (NumBits(avgjump)-1)/2 + dbadd;
  distbits += (distbits%2);
  // Set tau depending on whichever is faster, rounding tau' up or down.
  // Calculate the predicted running time for each kangaroo for rounding up and down.
  RR timepred1 = (2.0*sqrt(alpha*to_RR(U))/sqrt(to_RR(2.0*tau-1.0)) + to_RR(power2_ZZ(distbits))/tauprime )*(1.0+(tauprime - 1.0)/n);

  tau++;
  tauprime = (double)tau;  
  avgjump = CeilToZZ(2.0*sqrt(alpha*to_RR(U)*(2.0*tauprime-1.0)));
  distbits = (NumBits(avgjump)-1)/2 + dbadd;
  distbits += (distbits%2);
  RR timepred2 = (2.0*sqrt(alpha*to_RR(U))/sqrt(to_RR(2.0*tau-1.0)) + to_RR(power2_ZZ(distbits))/tauprime )*(1.0+(tauprime - 1.0)/n);
  
  tau++;
  tauprime = (double)tau;  
  avgjump = CeilToZZ(2.0*sqrt(alpha*to_RR(U)*(2.0*tauprime-1.0)));
  distbits = (NumBits(avgjump)-1)/2 + dbadd;
  distbits += (distbits%2);
  RR timepred3 = (2.0*sqrt(alpha*to_RR(U))/sqrt(to_RR(2.0*tau-1.0)) + to_RR(power2_ZZ(distbits))/tauprime )*(1.0+(tauprime - 1.0)/n);
  

  // Adjust tau depending on what would be faster.
  if((timepred2 < timepred1) && (timepred2 < timepred3))
    tau--;
  else if(timepred1 < timepred3)
    tau-=2;
  else;
  */

  // Set the average jump distance, the jump distances, 
  // and the starting position for tame.
  avgjump = RoundToZZ(sqrt(alpha*to_RR(U)*(2.0*tauratio-1.0)) - tauratio + 1.0);

  if(genus%3 == 1)
    headwind = (genus+2)/3;
  else
    headwind = genus/3;

  jumpdistance[63]=0;
  for(i=0; i<63; i++){
    jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 1 + headwind);
    
    jumpsum += jumpdistance[i];
    while( (jumpsum < (i-1)*(avgjump)) || (jumpsum > (i+1)*(avgjump)) ){
      jumpsum -= jumpdistance[i];
      jumpdistance[i] = (RandomBnd(2*avgjump-(genus+1)) + genus + 1 + headwind);
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
  below(E, Zero, tame);
  
  for(i=0; i<64; i++){
    below(jumpdistance[i], Zero, jumps[i]);
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
      baby_step_0_r2(tame, temp);
      tame = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to 0 as possible.
      if(tame.d2 < 0){
	reduce_basis2(tame);
	while(tame.d2 <= 0){
	  D=tame;
	  baby_step_2_r2(D, tame);
	  count++;
	}
	reduce_basis(D);
	for(j=0; j<count-1; j++){
	  tame = D;
	  baby_step_0_r2(tame, D);
	}
	tame = D;
	count=0;
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
	  if(check_inf2(order))
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
      baby_step_r1(wild, temp);
      wild = temp;
      
      // Correct to make sure that the 2-component of distance is 
      // as close to 0 as possible.
      if(wild.d2 < 0){
	reduce_basis2(wild);
	while(wild.d2 <= 0){
	  D=wild;
	  baby_step_2_r2(D, wild);
	  count++;
	}
	reduce_basis(D);
	for(j=0; j<count-1; j++){
	  wild = D;
	  baby_step_0_r2(wild, D);
	}
	wild = D;
	count=0;
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
	  if(check_inf2(order))
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

    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(tame.d2 < 0){
      reduce_basis2(tame);
      while(tame.d2 <= 0){
	D=tame;
	baby_step_2_r2(D, tame);
	count++;
      }
      reduce_basis(D);
      for(j=0; j<count-1; j++){
	tame = D;
	baby_step_0_r2(tame, D);
      }
      tame = D;
      count=0;
    }
    if(wild.d2 < 0){
      reduce_basis2(wild);
      while(wild.d2 <= 0){
	D=wild;
	baby_step_2_r2(D, wild);
	count++;
      }
      reduce_basis(D);
      for(j=0; j<count-1; j++){
	wild = D;
	baby_step_0_r2(wild, D);
      }
      wild = D;
      count=0;
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
  
  ZZ one = to_ZZ(1);

  ZZ R = extract2(order, one, factors, exponents, manext);
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
  
  if(order == R){
    cout<<"R     = "<<R<<endl;
    cout<<"h_x   = "<<order/R<<endl;
  }
  else{
    cout<<R<<" | R "<<R<<endl;
    cout<<"h_x | "<<order/R<<endl;
  }
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
  
  expgsjumps = 2.0*(alphaactual*to_RR(U)/to_RR(avgjump) + to_RR(avgjump)/(2.0*tauratio-1.0)) + to_RR(power2_ZZ(distbits+1))/tauratio;
  
  cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  cout<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
  cout<<"Total baby step kangaroo jumps:  "<<bsjumptotal<<endl;
  cout<<"Total giant step kangaroo jumps: "<<gsjumptotal<<endl;
  cout<<"Total kangaroo jumps:            "<<bsjumptotal+gsjumptotal<<endl;
  cout<<"Total distinguished points:      "<<distpoints<<endl;
  cout<<endl;

  //for(i=0; i<=genus; i++)
  //  cout<<"Frequency of headwind strength "<<i<<": "<<distdiffs[i]<<endl;
  //cout<<endl;
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
  
  tempnum = 2*E;
  below(tempnum, A);
  C = A;

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
void bsgs2(){
  ZZ R, order;
  infrastructure_ideal A, B, iI, iGS, C, D;
  cubic_ideal I, GS;
  ZZ GSdist, Idist;
  time_t set, bstime, gstime, checktime;
  long place1, place2, k;
  int j, count = 0;
  ZZ_pX s, s1, s2;
  cubic_ideal T;
  ZZ i;
  ZZ tempnum, wind;
  long wind0, wind1, wind2;

  headwind2(wind0, wind1, wind2);
  wind = to_ZZ(wind2);

  L = to_ZZ(1000000);

  reduce_basis(A);

  cout<<"Making Baby Steps ... \n";
  set = time(NULL);  

  for(i=0; i<L; i++){
   
    baby_step_0_r2(A, B);
    A = B;

    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(A.d2 < 0){
      reduce_basis2(A);
      while(A.d2 <= 0){
	D=A;
	baby_step_2_r2(D, A);
	count++;
      }
      reduce_basis(D);
      for(j=0; j<count-1; j++){
	A = D;
	baby_step_0_r2(A, D);
      }
      A = D;
      count=0;
    }
  }
  
  bstime = time(NULL);

  cout<<"\n"<<"Computed "<<L<<" Baby Steps in "<<bstime-set<<" seconds. \nMaking Giant Steps ... \n";

  tempnum = A.d0;
  do {
    below(tempnum, wind, iGS);
    tempnum--;
  } while(iGS.d2 != wind);  
  C = A;
  GSdist = iGS.d0;
  inf_to_cubic(iGS, GS);
  
  /*  inverse(I, C);
  if(I.d2 > wind){
    reduce_basis1(I);
    while(I.d2 >= wind){
      C=I;
      baby_step_1_r2(C, I);
    }
    reduce_basis(C);
    I = C;
  }
  else if(I.d2 < wind){
    reduce_basis2(I);
    while(I.d2 <= wind){
      C = I;
      baby_step_2_r2(C, I);
      count++;
    }
    reduce_basis(C);
    for(j=0; j<count-1; j++){
      I = C;
      baby_step_0_r2(I, C);
    }
    I = C;
    count=0;
  }
  else;
  
  cout<<I.d0<<" "<<I.d1<<" "<<I.d2<<endl;
  */

  for (i=0; i<L; i++){    
    /*
    // Search the hash table for the inverse of A.
    inverse(I, A);
    // Correct so that the 2-component of distance is 
    // as close to 0 as possible.
    if(I.d2 > 0){
      reduce_basis1(I);
      while(I.d2 >= 0){
	C=I;
	baby_step_1_r2(C, I);
      }
      reduce_basis(C);
      I = C;
    }
    else if(I.d2 < 0){
      reduce_basis2(I);
      while(I.d2 <= 0){
	C = I;
	baby_step_2_r2(C, I);
	count++;
      }
      reduce_basis(C);
      for(j=0; j<count-1; j++){
	I = C;
	baby_step_0_r2(I, C);
      }
      I = C;
      count=0;
    }
    else;
    */
    /*
    B = infrastructure_ideal(A);
    A = B*I;

    if(A.d2 < 0){
      reduce_basis2(A);
      while(A.d2 <= 0){
	C = A;
	baby_step_2_r2(C, A);
	count++;
      }
      reduce_basis(C);
      for(j=0; j<count-1; j++){
	A = C;
	baby_step_0_r2(A, C);
      }
      A = C;
      count=0;
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
    */
    
    
    B = infrastructure_ideal(A);
    giant_step_r2(A, B, GS, GSdist, wind); 
    // A = B*GS;

    if(A.d2 < 0){
      reduce_basis2(A);
      while(A.d2 <= 0){
	C = A;
	baby_step_2_r2(C, A);
	count++;
      }
      reduce_basis(C);
      for(j=0; j<count-1; j++){
	A = C;
	baby_step_0_r2(A, C);
      }
      A = C;
      count=0;
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
  }
  cout<<endl;
  gstime = time(NULL); 
  cout<<"Computed "<<L<<" Giant Steps in "<<gstime-bstime<<" seconds."<<endl;
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
 
  //         A) Baby Step, Giant Step or
  //         B) Pollard's Kangaroo
  if(rank == 0){
    if(bsgs)
      bsgs0(a, l);
    if(kang)
      kangaroo0(a, l); 
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
      bsgs2();
    if(kang){
      bsgs2();
      // kangaroo2();
    }
  }
}
