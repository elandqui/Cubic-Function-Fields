/************************************************************/
/* kangaroo.cc                                              */
/*                                                          */ 
/* This is the slave program for parallelized kangaroos.    */
/* This function makes the kangaroo jumps and checks for    */
/* distinguished points, writing them to a file.            */
/* The function terminates when the master program writes   */
/* a file containing the solution.                          */
/*                                                          */ 
/* Written by Eric Landquist                                */
/************************************************************/


#include "kangaroo.h"

NTL_CLIENT

/***********************************************************/
/* Get the splitting information for a set of polynomials. */
/***********************************************************/

void getSplitting(){

  fstream sfile;
  char *filename;
  bool closed = true;
  int block = 2*roonum + (rootype=='t' ? 0 : 1);
  vec_pair_ZZ_pX_long factors;
  int i;
  ZZ begin, end;
  time_t set1, time1;
  
  begin = ((block*(q-1))/m)+1;
  if(block == 0)
    begin--;
  end = ((block+1)*(q-1))/m;

  for(i=0; i<21; i++) 
    factordegs[i]=0;
  
  // Get the factors of f(x).
  factors = berlekamp(inv(LeadCoeff(G))*G);
  for(i=0; i<factors.length(); i++){
    factordegs[deg(factors[i].a)-1]++;
  }    
  factors = berlekamp(inv(LeadCoeff(H))*H);
  for(i=0; i<factors.length(); i++){
    factordegs[deg(factors[i].a)-1]++;
  }    
    
  set1 = time(NULL);

  // Loop through a block of monic irreducible degree lambda polynomials.
  SvCache.SetLength(2*lambda);
  Sv1(lambda, 1, begin, end);

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
  sfile<<splitting[1]<<"\n"<<splitting[2]<<"\n"<<i1time;
  sfile.close();  
}


/*************************************/
/* Make the jumps. Unit rank 0 case. */
/*************************************/

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

    i = vmap(roo);
    temp = cubic_ideal(roo);
    roo = temp*jumps[i];
    distance += jumpdistance[i];
    jumptotal++;

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
	found = distpoint(roo, distpoints, bsjumptotal, gsjumptotal, Zero, One);//, false);
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
  ZZ Stau = RoundToZZ(to_RR(q)/tau5());
  long dist0, dist1, dist2;
  ZZ wind0, wind2; 
  ZZ tempd1;
  //bool checkLattice = true;
  //ZZ icn = One;

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
    while(rep(coeff(roo.d, 0)) >= Stau){
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
	//if(checkLattice){
	found = distpoint(roo, distpoints, bsjumptotal, gsjumptotal, Zero, One);//, checkLattice);
	/*
	  if(!checkLattice){
	    // If we need to begin the lattice portion of the kangaroo
	    // algorithm, then re-initialize all the variables.
	    getLatticeRooInfo(icn);
	    bsjumptotal = Zero, gsjumptotal = Zero;
	    distbits2 = distbits/2;
	    /*************************/
	    /* Initialize the jumps. */
	    /*************************/
	    /*
	    for(i=0; i<64; i++){
	      jumps[i].s  = to_ZZ_pX(jumpsz[6*i]);
	      jumps[i].s1 = to_ZZ_pX(jumpsz[6*i+1]);
	      jumps[i].s2 = to_ZZ_pX(jumpsz[6*i+2]);
	      jumps[i].u  = to_ZZ_pX(jumpsz[6*i+3]);
	      jumps[i].v  = to_ZZ_pX(jumpsz[6*i+4]);
	      jumps[i].w  = to_ZZ_pX(jumpsz[6*i+5]);
	      }*/
	    /****************************/
	    /* Initialize the kangaroo. */
	    /****************************/
	/*
	    // This is a tame kangaroo.
	    if(rootype=='t'){
	      // Set the tame kangaroo.
	      distance = start + offset*roonum;
	      below(distance, icn, roo);
	    }
	    // This is a wild kangaroo.
	    else{
	      distance = offset*roonum;
	      below(distance, Zero, roo);
	    }
	  }
	  
	}
	else{
	  found = distpoint(roo, distpoints, bsjumptotal, gsjumptotal, checkLattice); //silly bear you are crazy and i love you!
	}*/
      }
    }
    
    /******************/
    /* Make the jump. */
    /******************/

    i = vmap(roo);
    temp = infrastructure_ideal(roo);
    giant_step_r2(roo, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    //roo = temp*jumps[i];
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
    gsjumptotal++;
  }
}

/*************************************/
/* Make the jumps. Unit rank 2 case. */
/*************************************/

void kangaroo2(ZZ &a, ZZ &l){
  int i, j, found = 0;
  ZZ dptest, roocoeff; // For distintinguished point checking.
  infrastructure_ideal roo;     // The kangaroo, tame or wild.
  infrastructure_ideal temp, B;
  cubic_ideal jumps[64];
  int distbits2 = distbits/2;
  ZZ distance = Zero;
  //ZZ jumptotal = ZZ(); // The number of jumps for this kangaroo.
  long distpoints = 0; // The number of distinguished points found.
  ZZ bsjumptotal = Zero, gsjumptotal = Zero;
  fstream dpfile;
  double tauratio = tau5();
  ZZ Stau;
  long dist0, dist1, dist2;
  ZZ wind0, wind2; 
  ZZ tempd1;
  ZZ eclass;
  // Make sure the roo hops in the correct equivalence class.
  if(rootype == 't')
    eclass = a;
  else
    eclass = Zero;
  //bool checkLattice = true;
  //ZZ icn = One;

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
	do{
	  baby_step_0_r2(roo, temp);
	  roo = temp;
	  reduce_basis2(temp);
	  while(temp.d2 <= 0){
	    B = temp;
	    baby_step_2_r2(B, temp);
	  }
	  roo = B;
	  reduce_basis(roo);
	} while(roo.d2 < 0);
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
	//if(checkLattice){
	found = distpoint(roo, distpoints, bsjumptotal, gsjumptotal, eclass, l);//, checkLattice);
	/*
	  if(!checkLattice){
	    // If we need to begin the lattice portion of the kangaroo
	    // algorithm, then re-initialize all the variables.
	    getLatticeRooInfo(icn);
	    bsjumptotal = Zero, gsjumptotal = Zero;
	    distbits2 = distbits/2;
	    /*************************/
	    /* Initialize the jumps. */
	    /*************************/
	    /*
	    for(i=0; i<64; i++){
	      jumps[i].s  = to_ZZ_pX(jumpsz[6*i]);
	      jumps[i].s1 = to_ZZ_pX(jumpsz[6*i+1]);
	      jumps[i].s2 = to_ZZ_pX(jumpsz[6*i+2]);
	      jumps[i].u  = to_ZZ_pX(jumpsz[6*i+3]);
	      jumps[i].v  = to_ZZ_pX(jumpsz[6*i+4]);
	      jumps[i].w  = to_ZZ_pX(jumpsz[6*i+5]);
	      }*/
	    /****************************/
	    /* Initialize the kangaroo. */
	    /****************************/
	/*
	    // This is a tame kangaroo.
	    if(rootype=='t'){
	      // Set the tame kangaroo.
	      distance = start + offset*roonum;
	      below(distance, icn, roo);
	    }
	    // This is a wild kangaroo.
	    else{
	      distance = offset*roonum;
	      below(distance, Zero, roo);
	    }
	  }
	  
	}
	else{
	  found = distpoint(roo, distpoints, bsjumptotal, gsjumptotal, checkLattice); //silly bear you are crazy and i love you!
	}*/
      }
    }
    
    /******************/
    /* Make the jump. */
    /******************/

    i = vmap(roo);
    temp = infrastructure_ideal(roo);
    giant_step_r2(roo, temp, jumps[i], jumpdistance[i], jumpdistance1[i], wind2);
    //roo = temp*jumps[i];
    // Correct to make sure that the 2-component of distance is 
    // as close to 0 as possible.
    if(roo.d2 < 0){
      do{     
	baby_step_0_r2(roo, temp);
	roo = temp;
	reduce_basis2(temp);
	while(temp.d2 <= 0){ 
	  B = temp;
	  baby_step_2_r2(B, temp);
	} 
	roo = B; 
	reduce_basis(roo);
      } while(roo.d2 < 0);
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

/******************/
/* Hash functions */
/******************/

// A function mapping a kangaroo to {0, 1, ..., 63}

inline int vmap(cubic_ideal &A){
  return(trunc_long(rep(coeff(A.u, 0)), 6));
}

// A function mapping a kangaroo to {0, 1, ..., 63}

inline int vmap(infrastructure_ideal &A){
  return(trunc_long(rep(coeff(A.mu1, 0)), 6));
}

/****************************************************/
/* Process a distinguished point. Unit rank 0 case. */
/****************************************************/

// Write the distinguished point to a file.
// Return 0 if a solution has not been found, and 1 if a solution has been found.
int distpoint(cubic_ideal &roo, ZZ &distance, long &dps, ZZ &jumps){

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

    dps = 0;
    
    return(0);
  }
  roofile.close();

  // Open the distinguished point file.
  dpfile1.open(alldpname, ios::app);
  if( dpfile1.is_open() ){
    
    // Write the distinguished points and distance to the dpoints.roo file.
    // Structure of dpoints.roo:
    // roo.s [return] distance [return] rootype [return] roonum [return]
    dpfile1<<"\n"<<roo.s<<"\n"<<distance;//<<"\n"<<rootype<<"\n"<<roonum;
    
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
  dpfile2<<dps<<"\n"<<jumps<<"\n"<<roo.s<<"\n"<<roo.s1<<"\n"<<roo.s2<<"\n"<<roo.u<<"\n"<<roo.v<<"\n"<<roo.w<<"\n"<<distance<<endl;
  dpfile2.close();

  return(0);
}

/**********************************************************/
/* Process a distinguished point. Unit rank 1 and 2 case. */
/**********************************************************/

// Write the distinguished point to a file.
// Return 0 if a solution has not been found, and 1 if a solution has been found.
int distpoint(infrastructure_ideal &roo, long &dps, ZZ &bsjumps, ZZ &gsjumps, ZZ &eclass, ZZ &l){//, bool &checkLattice){
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
    infrastructure_ideal shiftjump, temp, B;
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
	do{
	  baby_step_0_r2(roo, temp); 
	  roo = temp;
	  reduce_basis2(temp);
	  while(temp.d2 <= 0){
	    B = temp;
	    baby_step_2_r2(B, temp);
	  }
	  roo = B; 
	  reduce_basis(roo);
	} while(roo.d2 < 0);
      }          
    }
    while(roo.d0%l != eclass){
      temp = roo;
      baby_step_0_r2(temp, roo);
    }
    gsjumps++;

    dps = 0;

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
/****************************************************************/
/* Process a distinguished point in the lattice part of rank 2. */
/****************************************************************/

// Write the distinguished point to a file.
// Return 0 if a solution has not been found, and 1 if a solution has been found.
/*
int distpointLat(infrastructure_ideal &roo, long dps, ZZ &bsjumps, ZZ &gsjumps){

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

  ss << "reassignL_" << rootype << roonum << ".roo";
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
    long wind0, wind1, wind2;
    headwind(wind0, wind1, wind2);
    ZZ wind = to_ZZ(wind2);
    below(shift, wind, shiftjump);
    
    temp = infrastructure_ideal(roo);
    roo = temp*shiftjump;
    
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

    gsjumps++;

    // Create the distinguished point filenames for this kangaroo.
    std::stringstream ss;
    lastdpname = (char *)malloc(17*sizeof(char));
    
    ss << "last_dp_" << rootype << roonum << ".roo";
    ss >> lastdpname;
    
    std::stringstream ss2;
    alldpname = (char *)malloc(17*sizeof(char));
    
    ss2 << "Lall_dp_" << rootype << roonum << ".roo";
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
  dpfile2<<dps<<"\n"<<gsjumps<<"\n"<<bsjumps<<"\n"<<roo.mu0<<"\n"<<roo.mu1<<"\n"<<roo.mu2<<"\n"<<roo.nu0<<"\n"<<roo.nu1<<"\n"<<roo.nu2<<"\n"<<roo.d<<"\n"<<roo.d0<<"\n"<<roo.d1<<"\n"<<roo.d2<<endl;
  dpfile2.close();

  return(0);
}
*/
// Reads initialization information from ini.roo.

void getRooInfo(){
  fstream roofile;
  bool closed = true;
  //char *token;
  int i, j=0;

  while(closed){
    roofile.open("ini.roo", ios::in);
    if( roofile.is_open() ) closed = false;
    else {
      roofile.close();
      sleep(1);
    }  
  }

  roofile >> start; // = to_ZZ(token);
  roofile >> offset; // = to_ZZ(token);
  roofile >> distbits; // = atoi(token);

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

  if(rank == 0){
    roofile >> gz[0];
    roofile >> gz[1];
    roofile >> gz[2];
    roofile >> gz[3];
    roofile >> gz[4];
    roofile >> gz[5];
  }
  for(i=0; i<64; i++){
    roofile >> jumpsz[6*i];
    roofile >> jumpsz[6*i+1];
    roofile >> jumpsz[6*i+2];
    roofile >> jumpsz[6*i+3];
    roofile >> jumpsz[6*i+4];
    roofile >> jumpsz[6*i+5];
  }
  roofile.close();
}

// Reads initialization information from latticeini.roo.
/*
void getLatticeRooInfo(ZZ &icn){
  fstream roofile;
  bool closed = true;
  //char *token;
  int i, j=0;

  while(closed){
    roofile.open("latticeini.roo", ios::in);
    if( roofile.is_open() ) closed = false;
    else {
      roofile.close();
      sleep(1);
    }  
  }

  roofile >> start; // = to_ZZ(token);
  roofile >> icn;   
  roofile >> offset; // = to_ZZ(token);
  roofile >> distbits; // = atoi(token);

  for(i=0; i<64; i++){
    roofile >> jumpdistance[i]; // = to_ZZ(token);
  }

  for(i=0; i<64; i++){
    roofile >> jumpsz[6*i];
    roofile >> jumpsz[6*i+1];
    roofile >> jumpsz[6*i+2];
    roofile >> jumpsz[6*i+3];
    roofile >> jumpsz[6*i+4];
    roofile >> jumpsz[6*i+5];
  }
  roofile.close();
}
*/
int main(int argc, char *argv[]){
  int i;
  bool closed = true;
  fstream gofile, sfile, dpfile;
  ZZ a, l; // h = a (mod l).
 
  if(argc != 3){
    cout<<"Usage: ./kangaroo t/w i \n\n";
    return 0;
  }

  /****************************/
  /* Initialize the kangaroo. */
  /****************************/

  if((strncmp(argv[1], "t", 1)==0) || (strncmp(argv[1], "T", 1)==0) )
    rootype = 't';
  else if((strncmp(argv[1], "w", 1)==0) || (strncmp(argv[1], "W", 1)==0) )
    rootype = 'w';
  else{
    cout<<"Kangaroo type unknown.\n\n";
    return(0);
  }

  roonum = atoi(argv[2]);

  // Wait for the getS.roo file before doing Phase 1 work.
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
  // Set the unit rank.
  if (deg(f)%3){
    genus = deg(G) + deg(H) - 1;
    rank = 0;
  }
  else {
    genus = deg(G) + deg(H) - 2;
    if(q%3 == 1)
      rank = 2;
    else if(isCube(LeadCoeff(f)))
      rank = 1;
    else
      rank = 0;
  }

  // Create the distinguished point filenames for this kangaroo.
  std::stringstream ss;
  lastdpname = (char *)malloc(17*sizeof(char));

  ss << "last_dp_" << rootype << roonum << ".roo";
  ss >> lastdpname;

  std::stringstream ss2;
  alldpname = (char *)malloc(16*sizeof(char));

  ss2 << "all_dp_" << rootype << roonum << ".roo";
  ss2 >> alldpname;
  
  // Check to see if we are picking up where we left off.
  dpfile.open(lastdpname, ios::in);
  if( dpfile.is_open() )  {
    dpfile.close();
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
      if(IsOne(l))
	kangaroo2();
      else
	kangaroo2(a, l);
    }
    
  }
  else{
    if((roonum < 0) || (roonum >= m/2)){
      cout<<"Bad kangaroo number 0 <= i < m/2.\n\n";
      return(0);
    }

    closed = true;

    // Get the prime splitting info for the assigned block.
    // If Phase 1 is quick, the master will take care of it.
    if((q%3 == 1) || (lambda%2 == 0))
      getSplitting();
    
    // Wait for the go signal before going.
    while(closed){
      gofile.open("go.roo", ios::in);
      if( gofile.is_open() ) closed = false;
      else {
        gofile.close();
        sleep(1);
      }  
    }
    gofile.close();
    
    /*********************************/
    /* Get information from ini.roo. */
    /*********************************/
    getRooInfo();

    /*****************************/
    /* Get the kangaroo hopping! */
    /*****************************/
    if(rank == 0)
      kangaroo0();
    else if(rank == 1){
      init(2*genus);
      kangaroo1();
    }
    else{
      init(2*genus);
      smallOrder0(a, l);
      if(IsOne(l))
	kangaroo2();
      else
	kangaroo2(a, l);
    }
  }

  return(1);

}
