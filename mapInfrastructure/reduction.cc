#include "reduction.h"

NTL_CLIENT

// xi_alpha = b*rho + c*omega
inline int xideg(ZZ_pX b, ZZ_pX c){
  return(deg(b*rho + c*omega) - shift);
}
// eta_alpha = b*rho - c*omega
inline int etadeg(ZZ_pX b, ZZ_pX c){
  return(deg(b*rho - c*omega) - shift);
}

// The degree of the discriminant of the ideal
inline int discdeg(ZZ_pX mu1,ZZ_pX mu2, ZZ_pX nu1, ZZ_pX nu2){
  return (2*(deg((mu1*rho + mu2*omega)*(nu1*rho - nu2*omega) - (mu1*rho - mu2*omega)*(nu1*rho + nu2*omega)) - 2*shift));
}

inline ZZ_pX divxi(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX ximu, xinu, q;
  ximu = mu0 + RightShift(mu1*rho + mu2*omega, shift);
  xinu = nu0 + RightShift(nu1*rho + nu2*omega, shift);
  div(q, ximu, xinu);
  return q;
}

inline ZZ_pX diveta(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX etamu, etanu, q;
  etamu = mu0 + RightShift(mu1*rho - mu2*omega, shift);
  etanu = nu0 + RightShift(nu1*rho - nu2*omega, shift);
  div(q, etamu, etanu);
  return q;
}

inline ZZ_pX trunczeta(ZZ_pX a, ZZ_pX b, ZZ_pX c){
  return(2*a - RightShift(b*rho + c*omega, shift));
}

inline ZZ_pX aut(ZZ_pX a, ZZ_pX b, ZZ_pX c){
  return(a + RightShift(u*b*rho + u*u*c*omega,shift));
}

void reduce(){
  int xm, xn, em, en; 
  ZZ_pX temp0, temp1, temp2, q, etamu, etanu;
  ZZ_p lc1, lc2;

  //Step 2 in 4.6
  xm = xideg(mu1, mu2);
  xn = xideg(nu1, nu2);

  if(xm < xn){
    temp0 = mu0; temp1 = mu1; temp2 = mu2;
    mu0 = nu0;  mu1 = nu1; mu2 = nu2;
    nu0 = -temp0; nu1 = -temp1; nu2 = -temp2;
  }
  if(xm = xn){
    em = etadeg(mu1, mu2); en = etadeg(nu1, nu2);
    if(em < en){
      temp0 = mu0; temp1 = mu1; temp2 = mu2;
      mu0 = nu0; mu1 = nu1; mu2 = nu2;
      nu0 = -temp0; nu1 = -temp1; nu2 = -temp2;
    }
  }

  // Step 3.
  em = etadeg(mu1, mu2); en = etadeg(nu1, nu2);
    if(em >= en){
      // Step 3.1

      while( deg(sqr(nu1*rho) - sqr(nu2*omega)) - 2*shift > (double)discdeg(mu1, mu2, nu1, nu2)/2  ){
	q = divxi(mu0, mu1, mu2, nu0, nu1, nu2);
	temp0 = mu0; temp1 = mu1; temp2 = mu2;
	mu0 = nu0; mu1 = nu1; mu2 = nu2;
	nu0 = q*mu0-temp0; nu1 = q*mu1-temp1; nu2 = q*mu2-temp2;
      }
      
      // Step 3.2
      q = divxi(mu0, mu1, mu2, nu0, nu1, nu2);
      temp0 = mu0; temp1 = mu1; temp2 = mu2;
      mu0 = nu0; mu1 = nu1; mu2 = nu2;
      nu0 = q*mu0-temp0; nu1 = q*mu1-temp1; nu2 = q*mu2-temp2;

      // Step 3.3
      etamu = mu1*rho - mu2*omega;
      etanu = nu1*rho - nu2*omega;
      
      em = deg(etamu); en = deg(etanu);
      if(em==en){
	lc1=ZZ_p(LeadCoeff(etamu));
	lc2=ZZ_p(LeadCoeff(etanu));
	lc1*=inv(lc2);
	
	mu0 -= lc1*nu0; mu1 -= lc1*nu1; mu2 -= lc1*nu2; 
      }
    }
    
    //Step 4.
    while(etadeg(mu1, mu2) > 0){
      q = diveta(nu0, nu1, nu2, mu0, mu1, mu2);
      temp0 = nu0; temp1 = nu1; temp2 = nu2;
      nu0 = mu0; nu1 = mu1; nu2 = mu2;
      mu0*=q; mu1*=q; mu2*=q;
      mu2-=temp2; mu1-=temp1; mu0-=temp0;
    }
    
    // Step 5.
    lc1 = inv(to_ZZ_p(2));
    mu0 -= lc1*trunczeta(mu0, mu1, mu2);
    nu0 -= lc1*trunczeta(nu0, nu1, nu2);

    // Normalize so mu and nu are monic.
    lc1 = ZZ_p(LeadCoeff(mu0 + RightShift(mu1*rho - mu2*omega, shift)));
    if(lc1 != 1){
      lc1 = inv(lc1);
      mu0*=lc1; mu1*=lc1; mu2*=lc1;
    }
    lc1 = ZZ_p(LeadCoeff(nu0 + RightShift(nu1*rho - nu2*omega, shift)));
    if(lc1 != 1){
      lc1 = inv(lc1);
      nu0*=lc1; nu1*=lc1; nu2*=lc1;
    }
}

// This function finds the first fundamental unit.

void unit1(){
  vec_ZZ_pX C, I; //C stores theta and I stores mu and nu.
  ZZ_pX t0=to_ZZ_p(1), t1, t2; // theta, elements of the preperiod and period.
  int clen=0, ilen=0;

  // Step 2 is reducing {1, rho, omega}.
  // Step 3.
  do{
    // Step 3.1
    clen+=3; ilen+=6;
    C.SetLength(clen);
    I.SetLength(ilen);
    C[clen-3] = t0; C[clen-2]=t1; C[clen-1]=t2;
    I[ilen-6] = mu0; I[ilen-5]=mu1; I[ilen-4]=mu2;
    I[ilen-3] = nu0; I[ilen-2]=nu1; I[ilen-1]=nu2;

    // Step 3.2
    if(deg(autnu = aut(nu0, nu1, nu2))==1)
    cout<<"The neighbor of 1 is {"<< nu0-ZZ_p(LeadCoeff(autnu))<<", "<<nu1<<", "<<nu2<<"}."<<endl;
  else
    cout<<"The neighbor of 1 is {"<<mu0<<", "<<mu1<<", "<<mu2<<"}."<<endl;


  }while(0);

}

int main(){

  int i;
  ZZ_pX autnu;

  // Set the characterisitc.
  ZZ_p::init(to_ZZ(p));
  
  // The coefficients of the basis elements.
  mu0 = ZZ_pX(0,0), mu1 = ZZ_pX(0,1), mu2 = ZZ_pX(0,0); 
  nu0 = ZZ_pX(0,0), nu1 = ZZ_pX(0,0), nu2 = ZZ_pX(0,1);
  

  // Initialize rho and omega.
  // The coefficients of rho and omega.
  int coefr[6] = {3,-1,1,-2,3,1};
  int coefo[7] = {1,3,1,-3,-2,-1,1};
  for(i = 0; i < 3+shift; i++)
    rho += ZZ_pX(i, coefr[i]);
  for(i = 0; i < 4+shift; i++)
    omega += ZZ_pX(i, coefo[i]);

  cout<<"RHO: "<<rho<<endl;
  cout<<"OMEGA: "<<omega<<endl;
  
  // Run the reduction.
  reduce();

  cout<<"mu: {" << mu0<<", "<<mu1<<", "<<mu2<<"}"<<endl;
  cout<<"nu: {" << nu0<<", "<<nu1<<", "<<nu2<<"}"<<endl;

  //Test for true reduction.

  if(xideg(mu1, mu2) > xideg(nu1, nu2)) {
    if(etadeg(mu1, mu2) < 1){
      if( 0< etadeg(nu1, nu2)){
	if(deg(trunczeta(mu0,mu1,mu2))<1){
	  if(deg(trunczeta(nu0,nu1,nu2))<1){
	    cout<<"The basis {1, mu, nu} is reduced."<<endl;
	  }}}}}

  

  unit1();
    
}
