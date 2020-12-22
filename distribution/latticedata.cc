/* 
   latticedata.cc

   This file contains functions to compute the divisor class number, h, 
   and regulator, R, of several purely cubic function fields of unit rank 2.
   We will also gather information about the fundamental units and the structure
   of the infrastructure. Included are routines for applying the Baby Step-Giant 
   Step algorithm of Shanks and also the Kangaroo method of Pollard.
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
   Compile latticedata via ./makeld in ../PointCount.
   
   Make sure that the paths to the NTL libraries are
   specified correctly in the makefiles.
   
   Known Problems: 
   ---------------
   Needs to be written.
   
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

#include "latticedata.h"

NTL_CLIENT

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

  set = time(NULL);
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

  //rootime = time(NULL);

  // Extract R_S as a factor of h_0.
  R = extract(order, One, factors, exponents, manext);
  //checktime = time(NULL);

  //RR alphaactual = to_RR(abs(order-E))/to_RR(U);
  
  //expgsjumps = 2.0*(2.0*alphaactual*to_RR(U)/to_RR(avgjump) + to_RR(avgjump)/(2.0*(2.0*tauratio-1.0))) + to_RR(power2_ZZ(distbits+1))/tauratio;
  
  //cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  //cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  //cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  //cout<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
  //cout<<"Total baby step kangaroo jumps:  "<<bsjumptotal<<endl;
  //cout<<"Total giant step kangaroo jumps: "<<gsjumptotal<<endl;
  //cout<<"Total kangaroo jumps:            "<<bsjumptotal+gsjumptotal<<endl;
  //cout<<"Total distinguished points:      "<<distpoints<<endl;
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
  cubic_ideal jumps[64];
  ZZ jumpdistance[64];
  ZZ jumpdistance1[64];
  ZZ jumpsum = ZZ();
  int j, count = 0;

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

  //cout<<" Done. Beginning the jumping."<<endl;
  //cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  //cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  //cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  //cout<<"Distinguished points expected:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
  //set = time(NULL);
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

  //rootime = time(NULL);

  // Extract R_S as a factor of h_0.
  ZZ R = extract2(order, One, factors, exponents, manext);
  
  // If R == order, then we know that h_x = 1 and we are done.
  // Otherwise, the regulator is some multiple of R. Determine which here.
  //if(R != order){
  icn = order/R;
  latticeSearch();
  //}
  //checktime = time(NULL);

  //cout<<"Phase 4: Completed in "<<checktime-rootime<<" seconds."<<endl;

  //RR alphaactual = to_RR(abs(order-E))/to_RR(U);
    
  //expgsjumps = (2.0*alphaactual*to_RR(U)/to_RR(avgjump) + 2.0*to_RR(avgjump)/(2.0*tauratio-1.0)) + to_RR(power2_ZZ(distbits));
  //cout<<"For this example:"<<endl;
  //cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  //cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  //cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  //cout<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  //cout<<endl;
  //cout<<"Total baby step kangaroo jumps:  "<<bsjumptotal<<endl;
  //cout<<"Total giant step kangaroo jumps: "<<gsjumptotal<<endl;
  //cout<<"Total kangaroo jumps:            "<<bsjumptotal+gsjumptotal<<endl;
  //cout<<"Total distinguished points:      "<<distpoints<<endl;
  //cout<<endl;
}

// Kangaroo method for unit rank 2.
// Input: a - an integer
//        l - an integer such that h = a (mod l)
// Output: Prints h, |h-E|, and |h-E|/U, along with timing data.
// A generalization of Alg. 6.2.11 of [L09].
void kangaroo2(ZZ &a, ZZ &l){
  int i, j, found = 0, count = 0;
  int distbits; // The number of zero bits to determine 
                // if an ideal is distinguished.
  int distbits2;
  int dbadd = 4;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  ZZ avgjump;  // The average of the jumps.
  ZZ randjump; // The upper bound on the random jump distance.
  infrastructure_ideal temp, tame, wild;
  cubic_ideal jumps[64];
  ZZ jumpdistance[64];
  ZZ jumpdistance1[64];
  ZZ jumpsum = ZZ();
  
  ZZ matchdist = ZZ();

  ZZ bsjumptotal = ZZ(); // The number of baby step  jumps for each kangaroo.
  ZZ gsjumptotal = ZZ(); // The number of giant step jumps for each kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  long hashsize;
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

  //cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  //cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  //cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  //cout<<"Distinguished points expected:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
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
	  temp=tame;
	  baby_step_2_r2(temp, tame);
	  count++;
	}
	reduce_basis(temp);
	for(j=0; j<count-1; j++){
	  tame = temp;
	  baby_step_0_r2(tame, temp);
	}
	tame = temp;
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
	  temp=wild;
	  baby_step_2_r2(temp, wild);
	  count++;
	}
	reduce_basis(temp);
	for(j=0; j<count-1; j++){
	  wild = temp;
	  baby_step_0_r2(wild, temp);
	}
	wild = temp;
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
	temp=tame;
	baby_step_2_r2(temp, tame);
      }
      reduce_basis(temp);
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
	temp=wild;
	baby_step_2_r2(temp, wild);
      }
      reduce_basis(temp);
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

  //rootime = time(NULL);

  //cout<<"Phase 3: Completed in "<<rootime-set<<" seconds."<<endl;
  
  //cout<<" |h-E|  = "<<abs(order-E)<<endl; 
  //cout<<"|h-E|/U = "<<to_RR(abs(order-E))/to_RR(U)<<"\n"<<endl;
  
  // Extract R_S as a factor of h_0.
  

  R = extract2(order, One, factors, exponents, manext);
  
  // If R == order, then we know that h_x = 1 and we are done.
  // Otherwise, the regulator is some multiple of R. Determine which here.
  //if(R != order){
  icn = order/R;
  latticeSearch();
  //}
  //checktime = time(NULL);

  //cout<<"Phase 4: Completed in "<<checktime-rootime<<" seconds."<<endl;

  //RR alphaactual = to_RR(abs(order-E))/to_RR(U);  
  //expgsjumps = (2.0*alphaactual*to_RR(U)/to_RR(avgjump) + 2.0*to_RR(avgjump)/(to_double(l)*(2.0*tauratio-1.0))) + to_RR(power2_ZZ(distbits));
  //cout<<"For this example:"<<endl;
  //cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ(to_double(l)*(tauratio-1.0)*expgsjumps)<<endl;
  //cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  //cout<<"Expected kangaroo jumps:         "<<RoundToZZ(to_double(l)*tauratio*expgsjumps)<<endl;
  //cout<<"Expected distinguished points:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  //cout<<endl;
  //cout<<"Total baby step kangaroo jumps:  "<<bsjumptotal<<endl;
  //cout<<"Total giant step kangaroo jumps: "<<gsjumptotal<<endl;
  //cout<<"Total kangaroo jumps:            "<<bsjumptotal+gsjumptotal<<endl;
  //cout<<"Total distinguished points:      "<<distpoints<<endl;
  //cout<<endl;
}


// A hash function for unit rank 1 and 2.
// Input: A - a kangaroo  
// Output: An element of {0, 1, ..., 63}

inline int vmap(infrastructure_ideal &A){
  return(rep(coeff(A.mu1,0))%64);
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
// for unit ranks 1 and 2, kangaroo method.
// Input: A - the kangaroo
//        size - the size of the hashtable.
inline long roohash(infrastructure_ideal &A, long size){
  return( ((rep(coeff(A.nu1, 0))<<hashbits) + rep(coeff(A.mu0, 1)))%size );
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
  ZZ tempnum;
  infrastructure_ideal A, B, C;
  vec_pair_ZZ_pX_long factors;
  
  tempnum = 2*E;
  below(tempnum, A);
  C = A;

  insert(A, baby1, baby2, baby_deg, hashsize);

  for(i=1; i<L; i++){

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
  infrastructure_ideal A, B, C, D;
  vec_pair_ZZ_pX_long factors;
  ZZ tempnum;
  long dist0, dist1, dist2;
  ZZ wind2; 

  headwind2(dist0, dist1, dist2);
  wind2 = to_ZZ(dist2);

  below(E, Zero, A);
  C = A;
  
  insert(A, baby1, baby2, baby_deg, hashsize);

  for(i=1; i<L; i++){
    baby_step_0_r2(A, B);
    A = B;
    insert(A, baby1, baby2, baby_deg, hashsize);
    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(A.d2 < 0){
      while(A.d2 <= 0){
	reduce_basis2(A);
	D=A;	
	baby_step_2_r2(D, A);
	reduce_basis(A);
	insert(A, baby1, baby2, baby_deg, hashsize);
      }
      reduce_basis(D);
      A = D;
      if(A.d2 < 0){
	do{
	  baby_step_0_r2(A, D);
	  insert(D, baby1, baby2, baby_deg, hashsize);
	  A = D;
	  while(D.d2 <= 0){
	    reduce_basis2(D);
	    B = D;
	    baby_step_2_r2(B, D);
	    reduce_basis(D);
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
  ZZ i;
  ZZ tempnum, wind;
  long wind0, wind1, wind2;

  headwind2(wind0, wind1, wind2);
  wind = to_ZZ(wind2);

  reduce_basis(A);

  insert(A, baby1, baby2, baby_deg, hashsize);

  for(i=1; i<bsNum; i++){
    baby_step_0_r2(A, B);
    A = B;
    insert(A, baby1, baby2, baby_deg, hashsize);
    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(A.d2 < 0){
      while(A.d2 <= 0){
	reduce_basis2(A);
	D=A;	
	baby_step_2_r2(D, A);
	reduce_basis(A);
	insert(A, baby1, baby2, baby_deg, hashsize);
      }
      reduce_basis(D);
      A = D;
      if(A.d2 < 0){
	do{
	  baby_step_0_r2(A, D);

	  insert(D, baby1, baby2, baby_deg, hashsize);
	  A = D;
	  while(D.d2 <= 0){
	    reduce_basis2(D);
	    B = D;
	    baby_step_2_r2(B, D);
	    reduce_basis(D);
	    insert(D, baby1, baby2, baby_deg, hashsize);
	  }
	  A = B;
	  reduce_basis(A);
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

  //set = time(NULL);
  order = makeBabySteps1(baby1, baby2, baby_deg, hashsize, GS);
  //bstime = time(NULL);

  if(!IsZero(order)){
    if(order < (E-(genus+2))){ 
    //  gstime = time(NULL);
    //  cout<<"\n"<<"The x-Regulator is R_x = "<<2*order<<".\n";
    //  cout<<"The S-Regulator is R^S = "<<order<<".\n\n";
    }
    else{ 
      //gstime = time(NULL);
      goto end;
    }

    return;
  }

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
  //gstime = time(NULL); 
    
  if(found == 1)
    order = (baby_deg[place1] - B.d0)/2;
  else{
    order = (baby_deg[place2] - I.d0)/2;
  }

 end:

  //cout<<"Phase 3: Completed in "<<gstime-set<<" seconds.\n"<<endl;
  
  // Extract R^S as a factor of h_0.
  int i;

  R = extract(order, GS.d0/2, factors, exponents, manext);

  //checktime = time(NULL);

  //cout<<"Phase 4: Completed in "<<checktime-gstime<<" seconds."<<endl;

  return;
}

// Baby Step-Giant Step algorithm for class number and regulator
// computation in unit rank 2 infrastructure.
// Input: The function field
// Output: Prints h, |h-E|, |h-E|/U, R, and h_x, along with 
// timing data and a factorization of h.
// A generalization of Alg. 6.2.6 of [L09].
void bsgs2(ZZ &a, ZZ &l){
  infrastructure_ideal A, AI, BI, B, iI, iGS, C, D;
  cubic_ideal I, GS;
  ZZ Id0, Id1, GSd0, GSd1;
  //time_t set, bstime, gstime, checktime;
  long place1, place2, k;
  int found = 0;
  ZZ_pX s, s1, s2;
  ZZ tempnum, wind;
  long wind0, wind1, wind2;
  int reps = 0;

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

  //set = time(NULL);
  order = makeBabySteps2(baby1, baby2, baby_deg, hashsize, iGS);
  //bstime = time(NULL);
  inf_to_cubic(iGS, GS);
  GSd0 = iGS.d0;
  GSd1 = iGS.d1;

  if(!IsZero(order)){
    if(order < (E-genus)){ 
      //gstime = time(NULL);
    }
    else{ 
      //gstime = time(NULL);
      goto end;
    }
    return;
  }

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
    place1 = search(A, baby1, baby2, baby_deg, hashsize);
    if(place1 >= 0){
    check1:
      order = baby_deg[place1] - A.d0;
      //cout<<order<<endl;
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

    place2 = search(AI, baby1, baby2, baby_deg, hashsize);
    if(place2 >= 0){
    check2:
      order = baby_deg[place2] - AI.d0;
      //cout<<order<<endl;
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
    if(reps){
      cout<<f<<" Bad giant step!"<<endl;
      return;
    }
    reps = 1;
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
  
  if(found == 1)
    order = baby_deg[place1] - B.d0;
  else{
    order = baby_deg[place2] - AI.d0;
  }

 end:

  //cout<<"Order OK? "<<check_inf2(order)<<endl;

  //cout<<"Phase 3: Completed in "<<gstime-set<<" seconds.\n"<<endl;
  
  // Extract a multiple of R as a factor of h_0.
  int i;

  R = extract2(order, GSd0, factors, exponents, manext);

  // If R == order, then we know that h_x = 1 and we are done.
  // Otherwise, the regulator is some multiple of R. Determine which here.
  // Need to add in functionality to run the Kangaroo method if 
  // there is insufficient memory.
  //if(R != order){
  icn = order/R;
  if(icn > 10000){
    cout<<f<<" Big potential ideal class number! "<<icn<<endl;
    return;
  }
  if(!latticeSearch()){
    place1 = place2 = -1;
    goto loop;
  }
      
  //}
  //checktime = time(NULL);

  //cout<<"Phase 4: Completed in "<<checktime-gstime<<" seconds."<<endl;

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

int latticeSearch(){
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
      latticeRoo();
    }
    
    // Initialize the hashtables.
    for(long i=0; i<hashsize; i++){
      baby_deg[i] = N1;
      clear(baby1[i]);
      clear(baby2[i]);
    }
    
    //cout<<"Making "<<bsNum<<" Baby Steps ... \n";
    //set = time(NULL);
    latticeBS(bsNum, baby1, baby2, baby_deg, hashsize, iGS);
    //bstime = time(NULL);

    inf_to_cubic(iGS, GS);
    GSd0 = iGS.d0;
    GSd1 = iGS.d1;
    wind = iGS.d2;
    //cout<<"\nCompleted baby steps. Making giant steps."<<endl;

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
   
    //cout<<endl;
    
    if(place1 >= 0){
    check1:
      unit2 = A.d0 - baby_deg[place1];
      processUnit(unit2, One);
      return 1;
    }
    else if(place2 >= 0){
      unit2 = C.d0 - baby_deg[place2];
      // Just to make sure that we get the 2-distance of the unit correct.
      // There's an extremely small chance that A.d2 != C.d2.
      below(unit2, C.d2, A);
      int l = 1;
      while(!IsOne(A.d)){
	C.d2 += l;
	below(unit2, C.d2, A);
	//cout<<A.d0<<" "<<A.d1<<" "<<A.d2<<endl;
	if(l>0) { l++;  l = -l; }
	else { l--; l = -l; }
      }
      
      //cout<<"Unit found at "<<A.d0<<"(oo_1 - oo_0) + "<<A.d2<<"(oo_1 - oo_2).\n\n";
      // Given the unit that we've found, we search the lattice
      // for the unit with minimal positive 2-distance in order
      // to nail down the regulator and ideal class number exactly.
      if(processUnit(unit2, C.d2))
	 return 1;
      else
	return 0;
    }
    else{
      if(repeat){
	ZZ tempnum = iGS.d0 - 2*genus; 
	do{
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
	//R *= icn;   icn = 1;    
      }
      else
	if(!latticeRoo())
	  return 0;
    }
  }
  else{
    if(!latticeRoo())
      return 0;
  }
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

int latticeRoo(){
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

  //cout<<"Initializing ...";
  //cout<<flush;

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
  expdps = to_long(RightShift(RoundToZZ(expgsjumps), distbits));

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

  //cout<<" Done. Beginning the jumping."<<endl;
  //cout<<"Exp. baby step kangaroo jumps:   "<<RoundToZZ((tauratio-1.0)*expgsjumps)<<endl;
  // cout<<"Exp. giant step kangaroo jumps:  "<<RoundToZZ(expgsjumps)<<endl;
  //cout<<"Expected kangaroo jumps:         "<<RoundToZZ(tauratio*expgsjumps)<<endl;
  //cout<<"Distinguished points expected:   "<<RightShift(RoundToZZ(expgsjumps), distbits)<<endl;
  
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
        
	//printf("\r");
	//printf("Distinguished points: %ld", distpoints);
	//cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(tame.d0-matchdist);
	  below(order, icn, C);
	  if(IsOne(C.d)){
	    //cout<<endl<<endl;
	    //cout<<"Unit found at "<<order<<"(oo_1 - oo_0) + "<<icn<<"(oo_1 - oo_2).\n\n";
	    processUnit(order, icn);
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
	
	//printf("\r");
	//printf("Distinguished points: %ld", distpoints);
	//cout<<flush;
        
	// WOO HOO! There's a match!
	if(found){
	  order = abs(wild.d0-matchdist);
	  below(order, icn, C);
	  if(IsOne(C.d)){
	    //cout<<endl<<endl;
	    //cout<<"Unit found at "<<order<<"(oo_1 - oo_0) + "<<icn<<"(oo_1 - oo_2).\n\n";
	    processUnit(order, icn);
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
    giant_step_r2(wild, temp, jumps[i], jumpdistance[i], jumpdistance1[i],  wind2);
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
//        R - potential regulator; the exponent of the Jacobian.
//        unit2 - the 0-distance of a found unit.
//        dist2 - the 2-distance of a found unit.

int processUnit(ZZ &unit2, ZZ &dist2){
  int i, j; 
  ZZ d2factor;
  infrastructure_ideal C;
  a0 = unit2;
  a2 = dist2;
  e10 = R;
  for(j=1; j<dist2; j++){
    if(((icn%j) == 0) && ((dist2%j) == 0)){
      for(i = 0; i < dist2/j; i++){
	if(IsZero((unit2+i*R)%(dist2/j))){
	  unit2 += i*R; unit2 /= (dist2/j); 
	  d2factor = to_ZZ(j); 
	  below(unit2, d2factor, C);
	  if(IsOne(C.d)){
	    //printUnits(R, Zero, C.d0, C.d2);
	    e20 = C.d0; e22 = C.d2;
	    R *= C.d2;   icn /= C.d2; 
	    orthogonalize();
	    return 1;
	  }
	  unit2 *= (dist2/j); unit2 -= i*R;
	}
      }
    }
  }
  if(!IsZero(icn%dist2))
    return 0;
  
  e20 = unit2; e22 = dist2;
  //printUnits(R, Zero, unit2, dist2);
  R*=dist2;  icn /= dist2;   

  orthogonalize();

  return 1;
}

// Given a system of fundamental units, finds a system as orthogonal as possible.

void orthogonalize(){
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

// getinput() reads in the data from the input file, initializes a number 
// of variables for the desired computation, and constructs a polynomial 
// to generate the function field if needed.
// Input: input - the file that we are reading data from.
int getinput(char *input){
  FILE *inputptr;
  char thisLine[128];
  char token[64], value[64];
  int rand=2;
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
      lambdaplus = lambda = atoi(value);
    else if(strncmp(token, "bsgs:", 5)==0)
      bsgs = atoi(value);
    else if(strncmp(token, "kang:", 5)==0)
      kang = atoi(value);
    else if(strncmp(token, "memory:", 7)==0)
      memory = atoi(value);
    else if(strncmp(token, "degreeG:", 8)==0)
      dG = atoi(value);
    else if(strncmp(token, "degreeH:", 8)==0)
      dH = atoi(value);
    else if(strncmp(token, "irreducible:", 12)==0)
      irred = atoi(value);
    else if(strncmp(token, "samples:", 8)==0)
      samples = atoi(value);
    else{}
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

// Get a random polynomial for testing.
// G and H are both monic. If deg(H) = 1, then H(x) = x.

void getPolynomial(){
  int i, df;
  vec_pair_ZZ_pX_long factors;

  // Get a random degree degG polynomial.

  G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);

  if(irred){
    while(!DetIrredTest(G))
      G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
  }

  // Make sure G is squarefree.
  while(!IsOne(GCD(G, diff(G)))){
    G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
  }

  // We set H(x) = x if deg(H) = 1 since that simplifies the arithmetic.
  if(dH == 0){
    H = ZZ_pX(0,1);
    f = G;
  }
  else if(dH == 1){
    SetX(H);
    f = LeftShift(G, 2);
    // Make sure x does not divide G.
    if(!IsOne(GCD(G,H))){
      while( !IsOne(GCD(G, diff(G))) && !IsOne(GCD(G,H)) ){
	G = ZZ_pX(dG, 1) + random_ZZ_pX(dG);
      }
    }
  }
  else{
    H = ZZ_pX(dH, 1) + random_ZZ_pX(dH);
    if(irred){
      do{
	while(!DetIrredTest(H))
	  H = ZZ_pX(dH, 1) + random_ZZ_pX(dH);
      }while(G == H);
    }
    // Make sure G and H are relatively prime and H is squarefree.
    while(!IsOne(GCD(G,H)) || !IsOne(GCD(H, diff(H)))  ){
      H = ZZ_pX(dH, 1) + random_ZZ_pX(dH);
    }
    f = G*sqr(H);
  }
  
  df = deg(f);

  if(df > 20){
    cout<<"Error: The degree of D(x) must be less than 21.\n\n";
    return;
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
  for(i=0; i<df; i++)
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
  return;
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

  if(rank == 0){
    cout<<"Not working with unit rank 0.\n"<<endl;
    return(0);
  }

  cout<<"h"<<"\t"<<"|h-E|"<<"\t"<<"|h-E|/U     "<<"\t"<<"E"<<"\t"<<"U"<<"\t"<<"R"<<"\t"<<"h_x"<<"\t"<<"u10"<<"\t"<<"u12"<<"\t"<<"u20"<<"\t"<<"u22"<<"\t"<<"e10"<<"\t"<<"e20"<<"\t"<<"e22"<<"\t"<<"found0"<<"\t"<<"found2"<<"\t"<<"polynomial"<<"\t"<<"factors of h"<<endl;

  for(int i=0; i<samples; i++){
    getPolynomial();
    //SetX(H);
    //G = ZZ_pX(4,1) + ZZ_pX(3, 431) + ZZ_pX(2,891) + ZZ_pX(1,146)+ZZ_pX(0,795);
    //f = G*H*H;
    //cout<<f<<"\t";
    //cout<<flush;
    // Step 1: Compute an approximation E of the class number h
    //         and an integer L such that |h - E| < L^2.
    //set = time(NULL);
    lambda = lambdaplus;
    approxh(0);

    //time1 = time(NULL);
    
    //cout<<"Phase 1: Completed in "<<time1-set<<" seconds.\n\n";

    // Step 2: Use extra information about h in (E-L^2, E+L^2),
    //         such as its distribution in the interval
    //         or h mod r for small primes r.
    
    smallOrder0(a, l);

    time2 = time(NULL);

    // Step 3: Find h in the interval (E-L^2, E+L^2) via
    //         A) Baby Step, Giant Step or
    //         B) Pollard's Kangaroo
    if(rank == 1){
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
    // Print the data.
    cout<<order<<"\t"<<abs(order-E)<<"\t"<<to_RR(abs(order-E))/to_RR(U)<<"\t"<<E<<"\t"<<U<<"\t"<<R<<"\t"<<icn<<"\t"<<u10<<"\t"<<u12<<"\t"<<u20<<"\t"<<u22<<"\t"<<e10<<"\t"<<e20<<"\t"<<e22<<"\t"<<a0<<"\t"<<a2<<"\t"<<f;

    for(int j=0; j < factors.length(); j++){
      cout<<"\t"<<factors[j]<<"\t"<<exponents[j];
    }
    cout<<endl;
  
    //time3 = time(NULL);
    //cout<<"Total running time: "<<time3-set<<" seconds.\n\n";
  }
}
