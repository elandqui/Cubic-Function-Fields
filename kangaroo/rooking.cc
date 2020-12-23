/************************************************************/
/* rooking.cc                                               */
/*                                                          */ 
/* This is the master program for parallelized kangaroos.   */
/* This program takes in information about a cubic function */
/* field, computes an approximation for its class number,   */
/* an error bound, and information for the tame and wild    */
/* kangaroo herds. It periodically opens the distinguished  */
/* point file and checks for a match between the paths of a */
/* tame kangaroo and a wild kangaroo.                       */
/*                                                          */ 
/* Written by Eric Landquist                                */
/************************************************************/

#include "rooking.h"

NTL_CLIENT
// Finds an approximation E of h and a bound L such that
// |h-E|<L^2.

void approxh(){
  int s1, s2;
  int n, v, divlim; // v|n
  int a, i, j, l;
  RR Epr;
  RR C, S, T; // Temp variables.
  RR Q = to_RR(q);
  RR A;
  RR psi3; //, psi2, psi1;
  ZZ SvSum; // Used for psi3.
  ofstream sfile;
  char *filename;
  fstream roofile;
  bool closed = true;
  long split1, split2;
  time_t set1, time1;
  int i1time;
  ZZ zero = ZZ(), qm1 = q-1;

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

  //if((q%3 == 2) && (lambda%2 == 0))
  //  lambda++;

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

  // This will most likely be the longest running loop 
  // before the kangaroo portion.
 
  sfile.open("getS.roo", ios::out);
  sfile<<q<<"\n"<<G<<"\n"<<H<<"\n"<<m<<"\n"<<lambda<<endl;
  sfile.close();

  // The master will get the degree < lambda polynomials
  // The slaves will get the degree = lambda polynomials.
  set1 = time(NULL);
  for(v=1; v<=lambda-1; v++){
    splitting[0] = splitting[1] = splitting[2] = splitting[3] = 0;
    
    if((q%3 == 1)||(v%2 == 0)) {
      Sv1(v, 1, zero, qm1);
      Sv1(v, 3, zero, qm1);
    } 
    else{
      // Note: This must be fixed for q=2(mod 3)
      // since S_v(a) takes on four values in this case:
      // S_v(1), S_v(2), S_v(3), and S_v(6).
      Sv2(v, 1);
      Sv2(v, 2);
    }
    Epr*=(power(power(Q,2*v)/(power(Q,2*v) + power(Q,v) + 1.0), splitting[1])*power(power(Q,2*v)/(power(Q,2*v) -2.0*power(Q,v) + 1.0), splitting[2])*power(power(Q,2*v)/(power(Q,2*v) - 1.0), splitting[3]));
  }
  // If q = 2 (mod 3) and lambda = 1 (mod 2), then Sv2(1) and Sv2(2) are easy.
  if((q%3 == 2) && (lambda%2 == 1)) {
    splitting[0] = splitting[1] = splitting[2] = splitting[3] = 0;
    Sv2(lambda, 1);
    Sv2(lambda, 2);
    Epr*=(pow(power(Q,2*v)/(power(Q,2*v) + power(Q,v) + 1.0),to_RR(splitting[1]))*pow(power(Q,2*v)/(power(Q,2*v) -2.0*power(Q,v) + 1.0),to_RR(splitting[2]))); 
  }
  time1 = time(NULL);
  
  phase1time = time1-set1;

  // Should probably sleep for a little bit here, but that's OK.

  // Gather up the degree lambda results from the slaves if they
  // collected splitting information.
  if((q%3 == 1)||(lambda%2 == 0)) {
    splitting[0] = splitting[1] = splitting[2] = splitting[3] = 0;
    for(i=0; i<m; i++){
      std::stringstream ss;
      filename = (char *)malloc(sizeof(char)*15);
      ss << "block" << i << ".roo";
      ss >> filename;
      
      while(closed){ 
	ifstream roofile;
	roofile.open(filename, ios::in);
	if(roofile.is_open()){
	  roofile>>split1;
	  roofile>>split2;
	  roofile>>i1time;
	  
	  roofile.close();
	  closed=false;
	}
	else{
	  roofile.close();
	  sleep(1);
	}
      }
      
      closed = true;
      
      splitting[1]+=split1;
      splitting[2]+=split2;
      phase1time+=i1time;
    }

    // Get the S-values needed to compute psi3.
    SvCache[2*(lambda-1)] = 2*splitting[2]-splitting[1];

    v = lambda;
    Sv1(lambda, 3, zero, qm1);

    Epr*=(power(power(Q,2*v)/(power(Q,2*v) + power(Q,v) + 1.0), splitting[1])*power(power(Q,2*v)/(power(Q,2*v) -2.0*power(Q,v) + 1.0), splitting[2])*power(power(Q,2*v)/(power(Q,2*v) - 1.0), splitting[3]));
  }

  E = RoundToZZ(Epr);

  /************************************************/
  /* Approximate the square root of the error, L. */
  /************************************************/
  T = sqrt(Q);
  RR Tinv = inv(T);

  i = lambda+1;
  n = i+1;
  
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

  U = CeilToZZ(to_RR(E)*expm1(psi3)+1.0/2.0);
  
  // I'm going with this guy to speed things up.
  L = CeilToZZ(sqrt(alpha*to_RR(U)));
}

/*********************************************/
/* Get the kangaroo jumps: unit rank 0 case. */
/*********************************************/

void getRooJumps0(ZZ &a, ZZ &l){
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
      i=0;
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
  distbits = (NumBits(avgjump)-1)/2 + 3;// (NumBits(avgjump)-1)/2;
  distbits += (distbits%2);

  //cout<<avgjump<<" "<<distbits<<endl;

  /****************************************/
  /* Print start and distbits to roo.ini. */
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
  ofstream gofile;
  gofile.open("go.roo",ios::out);
  /*gofile<<"              /\ /l            "<<endl;
  gofile<<"             ((.Y(!            "<<endl;
  gofile<<"              \ |/             "<<endl;
  gofile<<"              /  6~6,          "<<endl;
  gofile<<"              \ _    +-.       "<<endl;
  gofile<<"               \`-=--^-'       "<<endl;
  gofile<<"                \ \            "<<endl;
  gofile<<"               _/  \           "<<endl;
  gofile<<"              (  .  Y          "<<endl;
  gofile<<"             /\"\ `--^--v--.   "<<endl;
  gofile<<"            / _ `--\"T~\/~\/   "<<endl;
  gofile<<"           / \" ~\.  !         "<<endl;
  gofile<<"     _    Y      Y./'          "<<endl;
  gofile<<"    Y^|   |      |~~7          "<<endl;
  gofile<<"    | l   |     / ./'          "<<endl;
  gofile<<"    | `L  | Y .^/~T            "<<endl;
  gofile<<"    |  l  ! | |/| |            "<<endl;
  gofile<<"    | .`\/' | Y | !            "<<endl;
  gofile<<"    l  \"~   j l j_L______     "<<endl;
  gofile<<"     \,____{ __\"~ __ ,\_,\_   "<<endl;
  gofile<<endl<<"Art courtesy of "<<endl;
  gofile<<"http://the.sunnyspot.org/asciiart/gallery/kangaroo.html"<<endl;*/
  gofile<<"Boing! Boing! Boing!"<<endl;
  
  gofile.close();

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

  // Write the go file.
  ofstream gofile;
  gofile.open("go.roo",ios::out);               
  /*gofile<<"              /\ /l            "<<endl;
  gofile<<"             ((.Y(!            "<<endl;
  gofile<<"              \ |/             "<<endl;
  gofile<<"              /  6~6,          "<<endl;
  gofile<<"              \ _    +-.       "<<endl;
  gofile<<"               \`-=--^-'       "<<endl;
  gofile<<"                \ \            "<<endl;
  gofile<<"               _/  \           "<<endl;
  gofile<<"              (  .  Y          "<<endl;
  gofile<<"             /\"\ `--^--v--.   "<<endl;
  gofile<<"            / _ `--\"T~\/~\/   "<<endl;
  gofile<<"           / \" ~\.  !         "<<endl;
  gofile<<"     _    Y      Y./'          "<<endl;
  gofile<<"    Y^|   |      |~~7          "<<endl;
  gofile<<"    | l   |     / ./'          "<<endl;
  gofile<<"    | `L  | Y .^/~T            "<<endl;
  gofile<<"    |  l  ! | |/| |            "<<endl;
  gofile<<"    | .`\/' | Y | !            "<<endl;
  gofile<<"    l  \"~   j l j_L______     "<<endl;
  gofile<<"     \,____{ __\"~ __ ,\_,\_   "<<endl;
  gofile<<endl<<"Art courtesy of "<<endl;
  gofile<<"http://the.sunnyspot.org/asciiart/gallery/kangaroo.html"<<endl;*/
  gofile<<"Boing! Boing! Boing!"<<endl;
  gofile.close();
}

/*********************************************/
/* Get the kangaroo jumps: unit rank 2 case. */
/*********************************************/

void getRooJumps2(ZZ &a, ZZ &l){
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

  if(m > 2){
    offset = avgjump/((m/2)-1);
    offset -= (offset%l);
  }

  // Set the jumps. 
  for(i=0; i<64; i++){
    below(jumpdistance[i], wind2, temp);
    while((temp.d0%l != wind0) && (temp.d2 != wind2) ){
      jumpdistance[i] += l;
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

  // Write the go file.
  ofstream gofile;
  gofile.open("go.roo",ios::out);               
  /*  gofile<<"              /\ /l            "<<endl;
  gofile<<"             ((.Y(!            "<<endl;
  gofile<<"              \ |/             "<<endl;
  gofile<<"              /  6~6,          "<<endl;
  gofile<<"              \ _    +-.       "<<endl;
  gofile<<"               \`-=--^-'       "<<endl;
  gofile<<"                \ \            "<<endl;
  gofile<<"               _/  \           "<<endl;
  gofile<<"              (  .  Y          "<<endl;
  gofile<<"             /\"\ `--^--v--.   "<<endl;
  gofile<<"            / _ `--\"T~\/~\/   "<<endl;
  gofile<<"           / \" ~\.  !         "<<endl;
  gofile<<"     _    Y      Y./'          "<<endl;
  gofile<<"    Y^|   |      |~~7          "<<endl;
  gofile<<"    | l   |     / ./'          "<<endl;
  gofile<<"    | `L  | Y .^/~T            "<<endl;
  gofile<<"    |  l  ! | |/| |            "<<endl;
  gofile<<"    | .`\/' | Y | !            "<<endl;
  gofile<<"    l  \"~   j l j_L______     "<<endl;
  gofile<<"     \,____{ __\"~ __ ,\_,\_   "<<endl;
  gofile<<endl<<"Art courtesy of "<<endl;
  gofile<<"http://the.sunnyspot.org/asciiart/gallery/kangaroo.html"<<endl;*/
  gofile<<"Boing! Boing! Boing!"<<endl;
  gofile.close();
}

/****************************************/
/* Reads the distinguished point files. */
/****************************************/

int readDPfile(ZZ &l){
  ifstream dpfile, lastdpfile;
  long numdps = 0;
  bool closed = true;
  long i=0, j;
  int matchloc;
  char *filename, *lastdpfilename;
  ZZ d1, d2; // dummy variables

  dp.SetLength(size);
  dists.SetLength(size);

  type = (char*)malloc(sizeof(char)*size);
  roonum = (int*)malloc(sizeof(int)*size);

  //cout<<mt<<" "<<mw<<endl;

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
    //cout<<filename<<endl;
    closed = true;
    while(closed){
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
	  //cout<<j<<" of "<<numdps<<endl;
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
	dpfile.close();
	closed = false;
      }
      else  {
	dpfile.close();
	sleep(1);
      }  
    }
  }

  //cout<<"Sorting!"<<endl;

  // Sort the distinguished point arrays.
  sort(0, i-1);

  //cout<<"Searching!"<<endl;

  // Search for matches.
  return search(i, matchloc, l);

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
      swap4 = roonum[i]; roonum[i] = roonum[j]; roonum[j] = swap4;
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

int search(int len, int &matchloc, ZZ &l){
  int i;

  for(i=1; i<len; i++){
    if(dp[i-1]==dp[i]){
      // Two kangaroos of the same type collided. Reassign one
      if((type[i-1]==type[i]) && (dists[i-1] == dists[i])){
        //cout<<i<<" "<<roonum[i-1]<<" "<<roonum[i]<<endl;
        // Reassign the higher-numbered kangaroo.
	if(roonum[i-1] > roonum[i])
	  reassign(type[i-1], roonum[i-1], l);
	else
	  reassign(type[i], roonum[i], l);
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

void reassign(char tw, int num, ZZ &l){
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

// latticeRoo() applies the Kangaroo method in the Jacobian lattice
// with 2-distance of icn. When a unit (the identity) is found, we can 
// quickly find the actual ideal class number and regulator.
// Input: icn - potential ideal class number
//            - the 2-distance of a known unit
//        R - potential regulator
//          - the exponent of the Jacobian
// Output: The actual ideal class number, icn, and regulator, R, 
// This function is not complete.
/*
void latticeRoo(ZZ &icn, ZZ &R, ZZ &e10, ZZ &e12, ZZ &e20, ZZ &e22){
*/
 
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
    else if(strncmp(token, "degreeG:", 8)==0)
      dG = atoi(value);
    else if(strncmp(token, "degreeH:", 8)==0)
      dH = atoi(value);
    else if(strncmp(token, "kangaroos:", 10)==0){
      m = atoi(value);

      if(m%2){
        cout<<"The number of kangaroos must be even.\n\n";
        return(0);
      }
      mt = mw = m/2;
    }
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
      factordegs[deg(factors[i].a)-1]++;
    } 
  }
  // Set the genus and unit rank.
  if ((dG+2*dH)%3) {
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
  if(rank == 1)
    tauratio = tau3();
  if(rank == 2)
    tauratio = tau5();

  fclose(inputptr);
  return(1);
}


int main(int argc, char *argv[]){
  int found=0;
  int waittime;
  ZZ a, l; // h = a (mod l)
  bool closed = true;
  fstream roofile, dpfile;
  int i;
  int dptotal = 0, iroodp;
  ZZ gsjumptotal = ZZ(), bsjumptotal = ZZ(), iroogsjump, iroobsjump;
  char *filename;
  time_t setk, timek, timeR;
  char *token;
  vec_ZZ factors;
  vec_long exponents;
  RR alphaactual;
  ZZ icn;//, e10, e12 = Zero, e20, e22 = One;

  if(argc != 2){
    cout<<"Usage: ./rooking inputfile \n\n";
    return 0;
  }
  else {
    // Read in the data file.
    if(!getinput(argv[1])) {
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
       return(0);
    }
  }

  // Check if this is picking up mid-stream. We'll lose timings if this happens,
  // but there's no good way around that.
  fstream getfile2;
  getfile2.open("ini.roo", ios::in);
  if(getfile2.is_open() ){
    ZZ avgjump = to_ZZ(0), ijump;
    smallOrder0(a, l);
    fstream getfile;
    getfile.open("getS.roo", ios::in);
    getfile >> q;
    getfile >> G;
    getfile >> H;
    getfile >> m;
    getfile >> lambda;
    getfile.close();
    f = G*sqr(H);
    mt = mw = m/2;
    
    getfile2 >> E;
    getfile2 >> L; // Junk variable at this point.
    getfile2 >> distbits;
    for(i=0; i<64; i++){
      getfile2 >> ijump;
      avgjump += ijump;
      if(rank == 2)
	getfile2 >> ijump;
    }
    avgjump /= 64;
    getfile2.close();
    if(rank == 0){
      L = RoundToZZ(to_RR(2*avgjump)/(sqrt(to_RR(l))*to_RR(m)));
      U = sqr(L);
    }
    else if (rank == 1){
      tauratio = tau3();
      alpha = alphahat();
      L = RoundToZZ(to_RR(2*avgjump)/(sqrt(to_RR(2*tauratio-1)*to_RR(m))));
      U = sqr(L);
    }
    else{
      tauratio = tau5();
      alpha = alphahat();
      U = RoundToZZ(4.0*sqr(to_RR(avgjump) + to_RR(l)*(tauratio - 1))/(to_RR(alpha)*to_RR(m)*to_RR(m)*to_RR(l)*to_RR(2*tauratio-1)));
      L = RoundToZZ(sqrt(to_RR(U)));
    }
    if(rank > 0)
      init(2*genus);
    getfile2.close();
    goto collect;
  }
  getfile2.close();

  // Step 1: Compute an approximation E of the class number h
  //         and an integer L such that |h - E| < L^2.
  approxh();
  

  // Find elements of small order then get the initialization info.
  smallOrder0(a, l);
  if(rank == 0){  
    getRooJumps0(a, l);
  }
  else if(rank ==1){
    init(2*genus);
    getRooJumps1();
  }
  else{
    init(2*genus);
    getRooJumps2(a, l);
  }

collect:

  // Compute the delay between checking the distinguished point file.
  // This will be calculated so that we expect roughly 10 checks on 
  // average, or once a day, whichever is less.
  // This assumes 1000 jumps are performed each second for r=0
  // and 
  // Better timing data should be taken to check this.
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

  cout<<"Wait Time = "<<waittime<<endl;  

  // Set the expected size of the dp array.
  size = to_long(RightShift(to_ZZ(waittime)*to_ZZ(m)*10*rate, distbits));

  setk = time(NULL);
  
  while(!found){
    // Read in the distinguished point file every waittime seconds.
    sleep(waittime);
    found = readDPfile(l);
    //cout<<"Found a match? "<<found<<endl;
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

  // A solution has been found! Write the solution to a file.
  while(closed){
    roofile.open("solution.txt", ios::out);
    if( roofile.is_open() ) closed = false;
    else{
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
  /*
  if(rank == 2) {
    // Print the first unit.
    if(IsZero(e12))
      roofile<<"Units: div(unit_1) = "<<e10<<"(oo_1 - oo_0)"<<endl;
    else
      roofile<<"Units: div(unit_1) = "<<e10<<"(oo_1 - oo_0) + "<<e12<<"(oo_1 - oo_2)"<<endl;
    // Print the second unit.
    if(IsOne(e22))
      roofile<<"       div(unit_2) = "<<e20<<"(oo_1 - oo_0) + (oo_1 - oo_2)"<<endl;
    else
      roofile<<"       div(unit_2) = "<<e20<<"(oo_1 - oo_0) + "<<e22<<"(oo_1 - oo_2)"<<endl;
    cout<<endl;
  }
  */
  /****************************************************************/
  /* Round up data from the individual distinguished point files. */
  /****************************************************************/
  token = (char *)malloc(sizeof(char)*100);
  for(i=0; i<mt; i++){
    std::stringstream ss;
    filename = (char *)malloc(sizeof(char)*18);
    ss << "last_dp_t" << i << ".roo";
    ss >> filename;
    // Open the distinguished point file.
    dpfile.open(filename,ios::in);
    if( dpfile.is_open() )  {
      dpfile >> token; iroodp = atoi(token); 
      dpfile >> token; iroogsjump = to_ZZ(token);
      if(rank > 0){
	dpfile >> token; iroobsjump = to_ZZ(token);
      }
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
    // Open the distinguished point file.
    dpfile.open(filename,ios::in);
    if( dpfile.is_open() )  {
      dpfile>>token; iroodp = atoi(token); 
      dpfile>>token; iroogsjump = to_ZZ(token); 
      if(rank > 0){
	dpfile >> token; iroobsjump = to_ZZ(token);
      }
      dpfile.close();

      dptotal += iroodp;
      gsjumptotal += iroogsjump;
      if(rank > 0)
	bsjumptotal += iroobsjump;
    }
    else dpfile.close();
  }
  /**********************************/
  /* Print the final data and quit. */
  /**********************************/
  roofile<<"Distinguished Points found:      "<<dptotal<<endl;
  roofile<<"Total kangaroo giant step jumps: "<<gsjumptotal<<endl;
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
}
