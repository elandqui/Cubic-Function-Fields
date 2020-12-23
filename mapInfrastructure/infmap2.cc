#include "infmap2.h"

NTL_CLIENT

// Compute an approximation of A^(1/3) with prec precision.
ZZ_pX cuberoot(ZZ_pX A, int prec){
  ZZ_pX C;         // the cube root
  ZZ_pX T = A;     // temp polynomial
  int i, j, r, s, t;
  ZZ_p temp, denom;
  int k = deg(A);
  vec_ZZ_p Dseq; 
  vec_ZZ_p cube;
  int cubedeg = k/3;
  int ind = cubedeg + prec+1;

  cube.SetLength(ind);
  cube[ind-1] = 1;

  Dseq.SetLength(k+1);

  for(i=0; i<=k; i++){
    Dseq[k-i] = coeff(T, k-i);
    if( !IsZero(Dseq[k-i]) )
      T-=ZZ_pX(deg(T), LeadCoeff(T));
    if( IsZero(coeff(A, 0)) )
      Dseq[0] = 0;
  }

  for(i=1; i<ind; i++){
    // First term
    if(k-i >= 0)
      temp = Dseq[k-i];
    else
      temp=0;

    //cout<<"temp0 = "<<temp<<endl;

    // Second term
    for(t=ind-i+3; t<=ind; t++){
      for(s=ind-i+2; s<t; s++){
	for(r=ind-i+1; r<s; r++){
	  if(((r+s+t)-(3*prec+3)) == (k-i))
	    temp -= 6*cube[r-1]*cube[s-1]*cube[t-1];
	}
      }
    }

    //cout<<"temp  = "<<temp<<endl;

    // Third term
    for(r=ind-i+1; r<=ind; r++){
      for(s=ind-i+1; s<=ind; s++){
	if(((2*r+s)-3*(prec+1) == (k-i)) && (r!=s))
	  temp -= 3*sqr(cube[r-1])*cube[s-1];
      }
    }

    //cout<<"temp2 = "<<temp<<endl;

    // Fourth term
    if(!(i%3))
      temp -= power(cube[cubedeg+prec-(i/3)], 3);

    //cout<<"temp3 = "<<temp<<endl;

    //Finally
    if(i==1){
      denom = 3*sqr(cube[ind-1]);
      for(j=1; j<P; j++){
	if(IsOne(denom*j))
	  denom=to_ZZ_p(j);
      }
    }
    temp*=denom;

    //cout<<"temp4 = "<<temp<<endl;

    cube[ind-i-1] = temp;
  }
  for(j=0; j<ind; j++)
    C+=ZZ_pX(j, cube[j]);

  return C;
}

void basis(int prec){

  rho = cuberoot(LeftShift(G, 2), prec);
  omega = cuberoot(LeftShift(sqr(G),1),prec);

  //cout<<"RHO   = "<<rho<<endl;
  //cout<<"OMEGA = "<<omega<<endl;

}

// zeta_alpha = 2a - b*rho - c*omega
inline int zetadeg(ZZ_pX a, ZZ_pX b, ZZ_pX c){
  return(deg( LeftShift(2*a, shift) - b*rho - c*omega) - shift-ddeg);
}

// xi_alpha = b*rho + c*omega
inline int xideg(ZZ_pX b, ZZ_pX c){
  //cout<<"xi: "<<b*rho + c*omega<<endl;
  return(deg(b*rho + c*omega) - shift-ddeg);
}
// eta_alpha = b*rho - c*omega
inline int etadeg(ZZ_pX b, ZZ_pX c){
  //cout<<"eta: "<<b*rho - c*omega<<endl;
  return(deg(b*rho - c*omega) - shift-ddeg);
}
// xi_alpha = b*rho + c*omega
inline int xideg1(ZZ_pX b, ZZ_pX c){
  ZZ_p U = to_ZZ_p(u);
  //cout<<"xi: "<<b*rho + c*omega<<endl;
  return(deg(b*U*rho + c*U*U*omega) - shift-ddeg);
}
// eta_alpha = b*rho - c*omega
inline int etadeg1(ZZ_pX b, ZZ_pX c){
  ZZ_p U = to_ZZ_p(u);
  //cout<<"eta: "<<b*rho - c*omega<<endl;
  return(deg(b*U*rho - c*U*U*omega) - shift-ddeg);
}
// xi_alpha = b*rho + c*omega
inline int xideg2(ZZ_pX b, ZZ_pX c){
  ZZ_p U = to_ZZ_p(u);
  //cout<<"xi: "<<b*rho + c*omega<<endl;
  return(deg(b*U*U*rho + c*U*omega) - shift-ddeg);
}
// eta_alpha = b*rho - c*omega
inline int etadeg2(ZZ_pX b, ZZ_pX c){
  ZZ_p U = to_ZZ_p(u);
  //cout<<"eta: "<<b*rho - c*omega<<endl;
  return(deg(b*U*U*rho - c*U*omega) - shift-ddeg);
}
// The degree of the discriminant of the ideal
inline int discdeg(ZZ_pX mu1,ZZ_pX mu2, ZZ_pX nu1, ZZ_pX nu2){
  return (2*(deg((mu1*rho + mu2*omega)*(nu1*rho - nu2*omega) - (mu1*rho - mu2*omega)*(nu1*rho + nu2*omega)) -2*shift));
}

inline int eltdeg(ZZ_pX a, ZZ_pX b, ZZ_pX c){
  //cout<<"("<<LeftShift(a, shift) + b*rho + c*omega<<")/"<<d<<endl;

  return(deg(LeftShift(a, shift) + b*rho + c*omega) - shift-ddeg);
}

inline ZZ_pX divxi(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX ximu, xinu, q;
  ximu = mu1*rho + mu2*omega;
  xinu = nu1*rho + nu2*omega;

  div(q, ximu, xinu);
  return q;
}
inline ZZ_pX divxi1(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX ximu, xinu, q;
  ZZ_p U = to_ZZ_p(u);
  ximu = mu1*U*rho + mu2*U*U*omega;
  xinu = nu1*U*rho + nu2*U*U*omega;

  div(q, ximu, xinu);
  return q;
}
inline ZZ_pX divxi2(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX ximu, xinu, q;
  ZZ_p U = to_ZZ_p(u);
  ximu = mu1*U*U*rho + mu2*U*omega;
  xinu = nu1*U*U*rho + nu2*U*omega;

  div(q, ximu, xinu);
  return q;
}
inline ZZ_pX diveta(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX etamu, etanu, q;

  etamu = mu1*rho - mu2*omega;
  etanu = nu1*rho - nu2*omega;
  
  div(q, etamu, etanu);

  return q;
}
inline ZZ_pX diveta1(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX etamu, etanu, q;
  ZZ_p U = to_ZZ_p(u);
  etamu = mu1*U*rho - mu2*U*U*omega;
  etanu = nu1*U*rho - nu2*U*U*omega;
  
  //cout<<etamu<<etanu<<endl;

  div(q, etamu, etanu);

  return q;
}
inline ZZ_pX diveta2(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX etamu, etanu, q;
  ZZ_p U = to_ZZ_p(u);
  etamu = mu1*U*U*rho - mu2*U*omega;
  etanu = nu1*U*U*rho - nu2*U*omega;
  
  //cout<<etamu<<etanu<<endl;

  div(q, etamu, etanu);

  return q;
}
inline ZZ_pX trunczeta(ZZ_pX a, ZZ_pX b, ZZ_pX c){
  //ZZ_pX q;
  
  return (d*RightShift((LeftShift(2*a, shift)-b*rho-c*omega)/d, shift));
  
//  div(q, LeftShift(2*a, shift) - b*rho - c*omega, LeftShift(d, shift));

  //return(q*d);
}
inline ZZ_pX trunczeta1(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U){
  ZZ_pX q;
  div(q, LeftShift(2*a, shift) - b*U*rho - c*U*U*omega, LeftShift(d, shift));
  return(q*d);
}
inline ZZ_pX trunczeta2(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U){
  ZZ_pX q;
  div(q, LeftShift(2*a, shift) - b*U*U*rho - c*U*omega, LeftShift(d, shift));
  return(q*d);
}
inline ZZ_pX aut(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U){
  return(a + RightShift(U*b*rho + U*U*c*omega,shift));
}

inline int autdeg(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U){
  return(deg(LeftShift(a,shift) + U*b*rho + U*U*c*omega) - shift - ddeg);
}

void reduce(ZZ_p half){
  int xm, xn, em, en; 
  ZZ_pX temp0, temp1, temp2, q, etamu, etanu;
  ZZ_p lc, lc1;
   
  //Step 2 in 4.6
  xm = xideg(mu1, mu2);
  xn = xideg(nu1, nu2);

  if(xm < xn){
    temp0 = mu0; temp1 = mu1; temp2 = mu2;
    mu0 = nu0;  mu1 = nu1; mu2 = nu2;
    nu0 = -temp0; nu1 = -temp1; nu2 = -temp2;
  }
  else if(xm == xn){
    em = etadeg(mu1, mu2); en = etadeg(nu1, nu2);
    if(em < en){
      temp0 = mu0; temp1 = mu1; temp2 = mu2;
      mu0 = nu0; mu1 = nu1; mu2 = nu2;
      nu0 = -temp0; nu1 = -temp1; nu2 = -temp2;
    }
  }
  else;

  // Step 3.
  em = etadeg(mu1, mu2); en = etadeg(nu1, nu2);
  if(em >= en){
    // Step 3.1
    
    while( deg(sqr(nu1*rho) - sqr(nu2*omega)) - 2*shift  > (double)discdeg(mu1, mu2, nu1, nu2)/2  ){
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
      lc=ZZ_p(LeadCoeff(etamu));
      lc1=ZZ_p(LeadCoeff(etanu));
      lc*=inv(lc1);
      
      mu0 -= lc*nu0; mu1 -= lc*nu1; mu2 -= lc*nu2; 
    }
  }   

  //Step 4.
  while(etadeg(nu1, nu2) < 0){
    q = divxi(mu0, mu1, mu2, nu0, nu1, nu2);
      
    temp0 = mu0; temp1 = mu1; temp2 = mu2;
    mu0 = nu0; mu1 = nu1; mu2 = nu2;
    nu0*=q; nu1*=q; nu2*=q;
    nu2-=temp2; nu1-=temp1; nu0-=temp0;
    //cout<<"Step 4a: "<<endl;
    //cout<<"mu: {" << mu0<<", "<<mu1<<", "<<mu2<<"}"<<endl;
    //cout<<"nu: {" << nu0<<", "<<nu1<<", "<<nu2<<"}"<<endl;
  }
  
  while(etadeg(mu1, mu2) >= 0){ 
    q = diveta(nu0, nu1, nu2, mu0, mu1, mu2);
    
    temp0 = nu0; temp1 = nu1; temp2 = nu2;
    nu0 = mu0; nu1 = mu1; nu2 = mu2;
    mu0*=q; mu1*=q; mu2*=q;
    mu2-=temp2; mu1-=temp1; mu0-=temp0;
  }  

  // Step 5.
  //cout<<zetadeg(mu0, mu1, mu2)<<" "<<zetadeg(nu0, nu1, nu2)<<endl;
  //cout<<mu0<<half*trunczeta(mu0, mu1, mu2)<<endl;
  //cout<<nu0<<half*trunczeta(nu0, nu1, nu2)<<endl;
  //cout<<(LeftShift(2*mu0, shift)-mu1*rho-mu2*omega)<<4*RightShift(d*((LeftShift(2*mu0, shift)-mu1*rho-mu2*omega)/d), shift)<<4*trunczeta(mu0, mu1, mu2, flag)<<d<<endl;
  //cout<<mu0<<mu1<<mu2<<nu0<<nu1<<nu2<<d<<endl;

  //cout<<"|zeta_mu|_0 = "<<zetadeg(mu0, mu1, mu2)<<", |zeta_nu|_0 = "<<zetadeg(nu0, nu1, nu2)<<endl;

  //cout<<"|zeta_mu|_0 = "<<zetadeg(mu0, mu1, mu2)<<", |zeta_nu|_0 = "<<zetadeg(nu0, nu1, nu2)<<endl;

  //if(zetadeg(mu0, mu1, mu2) >= 0)
  mu0 -= half*trunczeta(mu0, mu1, mu2);
  
  //if(zetadeg(nu0, nu1, nu2) >= 0)
  nu0 -= half*trunczeta(nu0, nu1, nu2);

  //cout<<"|zeta_mu|_0 = "<<zetadeg(mu0, mu1, mu2)<<", |zeta_nu|_0 = "<<zetadeg(nu0, nu1, nu2)<<", |nu| = "<<eltdeg(nu0, nu1, nu2)<<", |eta_nu| = "<<etadeg(nu1,nu2)<<endl;

//  en = etadeg(nu1,nu2);
//  if((eltdeg(nu0, nu1, nu2) == en) && (en == 0)){
//    nu0 -= d*LeadCoeff(LeftShift(nu0, shift) + nu1*rho + nu2*omega);
//  }

  // Normalize so mu and nu are monic.
  lc = ZZ_p(LeadCoeff(mu0 + RightShift(mu1*rho - mu2*omega, shift)));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    mu0*=lc; mu1*=lc; mu2*=lc;
  }
  
  lc = ZZ_p(LeadCoeff(nu0 + RightShift(nu1*rho - nu2*omega, shift)));
  if(!(IsOne(lc)||IsZero(lc))){
    lc = inv(lc);
    nu0*=lc; nu1*=lc; nu2*=lc;
  }
  
}

// Computes a 1-reduced basis of {1, mu, nu}.
void reduce1(ZZ_p half, ZZ_p U){
  int xs, xt, es, et; 
  ZZ_pX temp0, temp1, temp2, q, etasig, etatau;
  ZZ_p lc, lc1;
  double disc_deg;

  //Step 2 in 4.6
  xs = xideg1(sig1, sig2);
  xt = xideg1(tau1, tau2);

  if(xs < xt){
    temp0 = sig0; temp1 = sig1; temp2 = sig2;
    sig0 = tau0;  sig1 = tau1; sig2 = tau2;
    tau0 = -temp0; tau1 = -temp1; tau2 = -temp2;
  }
  if(xs == xt){
    es = etadeg1(sig1, sig2); et = etadeg1(tau1, tau2);
    if(es < et){
      temp0 = sig0; temp1 = sig1; temp2 = sig2;
      sig0 = tau0; sig1 = tau1; sig2 = tau2;
      tau0 = -temp0; tau1 = -temp1; tau2 = -temp2;
    }
  }

  // Step 3.
  es = etadeg1(sig1, sig2); et = etadeg1(tau1, tau2);
    if(es >= et){
      // Step 3.1
      //cout<<"nu "<<tau1<<tau2<<endl;
      disc_deg = (double)discdeg(sig1, sig2, tau1, tau2)/2;
       while( deg(LeftShift(sqr(tau1),1)*U*U*omega - sqr(tau2)*G*U*rho) - shift  >  disc_deg){
	q = divxi1(sig0, sig1, sig2, tau0, tau1, tau2);
	temp0 = sig0; temp1 = sig1; temp2 = sig2;
	sig0 = tau0; sig1 = tau1; sig2 = tau2;
	tau0 = q*sig0-temp0; tau1 = q*sig1-temp1; tau2 = q*sig2-temp2;  
      }

      // Step 3.2
      q = divxi1(sig0, sig1, sig2, tau0, tau1, tau2);
      temp0 = sig0; temp1 = sig1; temp2 = sig2;
      sig0 = tau0; sig1 = tau1; sig2 = tau2;
      tau0 = q*sig0-temp0; tau1 = q*sig1-temp1; tau2 = q*sig2-temp2;

      // Step 3.3
      etasig = sig1*U*rho - sig2*U*U*omega;
      etatau = tau1*U*rho - tau2*U*U*omega;
      
      es = deg(etasig); et = deg(etatau);
      if(es==et){
	lc=ZZ_p(LeadCoeff(etasig));
	lc1=ZZ_p(LeadCoeff(etatau));
	lc*=inv(lc1);
	
	sig0 -= lc1*tau0; sig1 -= lc1*tau1; sig2 -= lc1*tau2; 
      }
    }

    //cout<<"sig: {" <<sig0<<", "<<sig1<<", "<<sig2<<"}/"<<d<<endl;
    //cout<<"tau: {" <<tau0<<", "<<tau1<<", "<<tau2<<"}/"<<d<<endl;

    //Step 4.
    while(etadeg1(tau1, tau2) < 0){
      q = divxi1(sig0, sig1, sig2, tau0, tau1, tau2);
      
      temp0 = sig0; temp1 = sig1; temp2 = sig2;
      sig0 = tau0; sig1 = tau1; sig2 = tau2;
      tau0*=q; tau1*=q; tau2*=q;
      tau2-=temp2; tau1-=temp1; tau0-=temp0;
      //cout<<"Step 4a: "<<endl;
      //cout<<"mu: {" << mu0<<", "<<mu1<<", "<<mu2<<"}"<<endl;
      //cout<<"nu: {" << nu0<<", "<<nu1<<", "<<nu2<<"}"<<endl;
    }

    while(etadeg1(sig1, sig2) >= 0){
      q = diveta1(tau0, tau1, tau2, sig0, sig1, sig2);

      temp0 = tau0; temp1 = tau1; temp2 = tau2;
      tau0 = sig0; tau1 = sig1; tau2 = sig2;
      sig0*=q; sig1*=q; sig2*=q;
      sig2-=temp2; sig1-=temp1; sig0-=temp0;
    }
    // Step 5.

    sig0 -= half*trunczeta1(sig0, sig1, sig2, U);
    tau0 -= half*trunczeta1(tau0, tau1, tau2, U);

    // et = etadeg(tau1,tau2);
    //if((eltdeg(tau0, tau1, tau2) == et) && (et == 0)){
    //  tau0 -= d*LeadCoeff(LeftShift(tau0, shift) + tau1*rho + tau2*omega);
    // }

}

// Computes a 2-reduced basis of {1, mu, nu}.
void reduce2(ZZ_p half, ZZ_p U){
  int xs, xt, es, et; 
  ZZ_pX temp0, temp1, temp2, q, etasig, etatau;
  ZZ_p lc, lc1;

  //Step 2 in 4.6
  xs = xideg2(sig1, sig2);
  xt = xideg2(tau1, tau2);

  if(xs < xt){
    temp0 = sig0; temp1 = sig1; temp2 = sig2;
    sig0 = tau0;  sig1 = tau1; sig2 = tau2;
    tau0 = -temp0; tau1 = -temp1; tau2 = -temp2;
  }
  if(xs == xt){
    es = etadeg2(sig1, sig2); et = etadeg2(tau1, tau2);
    if(es < et){
      temp0 = sig0; temp1 = sig1; temp2 = sig2;
      sig0 = tau0; sig1 = tau1; sig2 = tau2;
      tau0 = -temp0; tau1 = -temp1; tau2 = -temp2;
    }
  }

  // Step 3.
  es = etadeg2(sig1, sig2); et = etadeg2(tau1, tau2);
    if(es >= et){
      // Step 3.1
      //cout<<"nu "<<tau1<<tau2<<endl;
       while( deg(LeftShift(sqr(tau1),1)*U*omega - sqr(tau2)*G*U*U*rho) - shift  > (double)discdeg(sig1, sig2, tau1, tau2)/2  ){
	q = divxi2(sig0, sig1, sig2, tau0, tau1, tau2);
	temp0 = sig0; temp1 = sig1; temp2 = sig2;
	sig0 = tau0; sig1 = tau1; sig2 = tau2;
	tau0 = q*sig0-temp0; tau1 = q*sig1-temp1; tau2 = q*sig2-temp2;  
      }

      // Step 3.2
      q = divxi2(sig0, sig1, sig2, tau0, tau1, tau2);
      temp0 = sig0; temp1 = sig1; temp2 = sig2;
      sig0 = tau0; sig1 = tau1; sig2 = tau2;
      tau0 = q*sig0-temp0; tau1 = q*sig1-temp1; tau2 = q*sig2-temp2;

      // Step 3.3
      etasig = sig1*U*U*rho - sig2*U*omega;
      etatau = tau1*U*U*rho - tau2*U*omega;
      
      es = deg(etasig); et = deg(etatau);
      if(es==et){
	lc=ZZ_p(LeadCoeff(etasig));
	lc1=ZZ_p(LeadCoeff(etatau));
	lc*=inv(lc1);
	
	sig0 -= lc1*tau0; sig1 -= lc1*tau1; sig2 -= lc1*tau2; 
      }
    }

    //cout<<"sig: {" <<sig0<<", "<<sig1<<", "<<sig2<<"}/"<<d<<endl;
    //cout<<"tau: {" <<tau0<<", "<<tau1<<", "<<tau2<<"}/"<<d<<endl;

    //Step 4.
    while(etadeg2(tau1, tau2) < 0){
      q = divxi2(sig0, sig1, sig2, tau0, tau1, tau2);
      
      temp0 = sig0; temp1 = sig1; temp2 = sig2;
      sig0 = tau0; sig1 = tau1; sig2 = tau2;
      tau0*=q; tau1*=q; tau2*=q;
      tau2-=temp2; tau1-=temp1; tau0-=temp0;
      //cout<<"Step 4a: "<<endl;
      //cout<<"mu: {" << mu0<<", "<<mu1<<", "<<mu2<<"}"<<endl;
      //cout<<"nu: {" << nu0<<", "<<nu1<<", "<<nu2<<"}"<<endl;
    }

    while(etadeg2(sig1, sig2) >= 0){
      q = diveta2(tau0, tau1, tau2, sig0, sig1, sig2);

      temp0 = tau0; temp1 = tau1; temp2 = tau2;
      tau0 = sig0; tau1 = sig1; tau2 = sig2;
      sig0*=q; sig1*=q; sig2*=q;
      sig2-=temp2; sig1-=temp1; sig0-=temp0;
    }
    // Step 5.

    sig0 -= half*trunczeta2(sig0, sig1, sig2, U);
    tau0 -= half*trunczeta2(tau0, tau1, tau2, U);

    // et = etadeg(tau1,tau2);
    //if((eltdeg(tau0, tau1, tau2) == et) && (et == 0)){
    //  tau0 -= d*LeadCoeff(LeftShift(tau0, shift) + tau1*rho + tau2*omega);
    // }

}

// This function will compute the triangle basis of an ideal.
// The input is an ideal A = {1, mu, nu}/d = [d, mu0, mu1, ..., nu2]
// The return value is the vector [s, s', u, s", v, w] where
// J = {1, s'(u + rho), s"(v + w*rho + omega)}/s}

void triangle(vec_ZZ_pX A, vec_ZZ_pX &J){
  ZZ_pX a1, b1, a, b, t, u, w;

  cout<<A[0]<<A[1]<<A[2]<<A[3]<<A[4]<<A[5]<<A[6]<<endl;

  J[0] = A[0];
  XGCD(J[3], a1, b1, A[3], A[6]);
  J[1] = (A[2]*A[6] - A[5]*A[3])/J[3];

  if( !IsZero(t = (a1*A[2] + b1*A[5])%J[3]) ) {
    if( !IsZero(a = J[1]%J[3]) )
      t /= a;
     else    
       t = ZZ_pX();  // t = 0.
  }
  else
    t = ZZ_pX();  // t = 0.

  a = a1 - t*A[6]/J[3];
  b = b1 + t*A[3]/J[3];
  u = J[2] = (A[1]*A[6] - A[4]*A[3])/(J[1]*J[3]);
  J[4] = (a*A[1] + b*A[4])/J[3];
  w = J[5] = (a*A[2] + b*A[5])/J[3];

  J[2] %= (A[0]/J[1]);
  J[5] %= J[1];
  J[4] += u*(J[5]-w);
  J[4] %= (A[0]/J[3]);

  // Normalize s, s', and s".
  J[0]*=inv(LeadCoeff(J[0]));
  J[1]*=inv(LeadCoeff(J[1]));
  J[3]*=inv(LeadCoeff(J[3]));

  cout<<J[0]<<J[1]<<J[2]<<J[3]<<J[4]<<J[5]<<endl;

  return;
}

// This function finds J = A*B.
// J is in a triangle basis.
// Format: {s, s'(u + rho), s"(v + w rho + omega)}
// [s, s', u, s", v, w]
// Returns the degree of the factor divided out.

int compose(vec_ZZ_pX A, vec_ZZ_pX B, vec_ZZ_pX &J){

  ZZ_pX g, temp;
  vec_ZZ_pX P;
  P.SetLength(6);
  int divdeg=0;

  // If (s_A, s_B) = 1, then we do the easiest composition.
  if(IsOne(GCD(A[0], B[0]))){
    return relprimecomp(A, B, J);
  }

  else{
    XGCD(g, P[0], P[1], A[3]*B[3]*(A[4]+B[4]+LeftShift(A[5]*B[5],1)), B[1]*A[3]*(B[2]+LeftShift(A[5], 1)) );
    if(!IsOne(g)){
      XGCD(g, P[2], temp, B[0]*A[3], g);
      P[0]*=temp;
      P[1]*=temp;
      if(!IsOne(g)){
	XGCD(g, P[3], temp, A[1]*B[3]*(A[2]+LeftShift(B[5], 1)), g);
	P[0]*=temp;
	P[1]*=temp;
	P[2]*=temp;
	if(!IsOne(g)){
	  XGCD(g, P[4], temp, LeftShift(A[1]*B[1], 1), g);
	  P[0]*=temp;
	  P[1]*=temp;
	  P[2]*=temp;
	  P[3]*=temp;
	  if(!IsOne(g)){
	    XGCD(g, P[5], temp, A[0]*B[3], g);
	    P[0]*=temp;
	    P[1]*=temp;
	    P[2]*=temp;
	    P[3]*=temp;
	    P[4]*=temp;
	  }
	}
      }
    }

    // cout<<g<<endl;

    if (IsZero(coeff(A[0], 0)) && IsZero(coeff(B[0], 0)))
      SetX(temp);
    else
      temp = ZZ_pX(0,1);

    if(IsOne(g))
      return primcomp(A, B, P, 0, J);
    else if(g == A[3]*B[3]*temp)
      return primcomp(A, B, P, 1, J);
    else {
      //cout<<"Dividing primes: g = "<<g<<endl;
      vec_ZZ_pX I1, I2;
      I1.SetLength(6);
      I2.SetLength(6);
      divdeg += (deg(A[0]) + deg(B[0]) - deg(g) - divprimes(A, B, g, I1));
      //cout<<"Divprimes returns A, B, and I.\n";
      divdeg += compose(A, B, I2);
      //cout<<"Compose returns A*B=I.\n";
      divdeg += compose(I1, I2, J);	
      divdeg -= deg(J[0]);
      return (divdeg);
    }
  }
}

// This function finds J = A*B, where (s_A, s_B) = 1.
// J is in a triangle basis.
// Format: {s, s'(u + rho), s"(v + w rho + omega)}
// [s, s', u, s", v, w]

int relprimecomp(vec_ZZ_pX A, vec_ZZ_pX B, vec_ZZ_pX &J){
  
  ZZ_pX g, ra, rb, sa, sb, u, w;

  J[0] = A[0]*B[0];  // s_3 = s_1 s_2
  J[1] = A[1]*B[1];  // s_3' = s_1' s_2'
  J[3] = A[3]*B[3];  // s_3" = s_1" s_2"
  
  // Run CRT on the u's, v's, and w's.
  XGCD(g, ra, rb, sa = A[0]/A[1], sb = B[0]/B[1]);
  u = J[2] = (A[2]*rb*sb + B[2]*ra*sa)%(sa*sb);

  XGCD(g, ra, rb, A[1], B[1]);
  w = J[5] = (A[5]*rb*B[1] + B[5]*ra*A[1])%(A[1]*B[1]);

  XGCD(g, ra, rb, sa = A[0]/A[3], sb = B[0]/B[3]);
  J[4] = ((A[4] + A[2]*(J[5]-A[5]))*rb*sb + (B[4] + B[2]*(J[5]-B[5]))*ra*sa)%(sa*sb);

  // "Reduce."

  J[2] %= (J[0]/J[1]);
  J[5] %= J[1];
  J[4] += u*(J[5]-w);
  J[4] %= (J[0]/J[3]);

  //cout<<J[0]<<J[1]<<J[2]<<J[3]<<J[4]<<J[5]<<endl;

  return 0;
}

// This function finds J such that J = A^2.
// The ideal A must have a triangle basis. 
// The ideal J has a triangle basis.
// Return the degree of s_Gs_H/s_G's".

int square(vec_ZZ_pX A, vec_ZZ_pX &J) {
 
  ZZ_pX sg, sh, sg1, sgsh, sonsgh, sonsgh2, umod, g, a, b, sa, sb, temp, temp2;
  ZZ_pX y, z, u, w;
  ZZ_pX f;
  int Fdeg;

  // cout<<A[0]<<A[1]<<A[2]<<A[3]<<A[4]<<A[5]<<endl;

  GCD(sg, A[0], G);
  //  GCD(sh, A[0], H);
  if(IsZero(ConstTerm(A[0])))
    SetX(sh);
  else
    sh = ZZ_pX(0,1);
  GCD(sg1, A[1], G);

  sgsh = sg*sh;
  sonsgh = A[0]/sgsh;
  sqr(sonsgh2, sonsgh);
  umod = sonsgh*sg1/A[1];
  
  f = sgsh/(sg1*A[3]);
  //cout<<F<<endl;

  XGCD(g, a, y, umod, 3*sqr(A[2]));
  temp2 = 2*A[4] + LeftShift(sqr(A[5]),1);
  XGCD(g, a, z, sonsgh, temp2);
  //cout<<g<<endl;
  while(!IsOne(g)){
    //f = ZZ_pX(0,1);
    A[5] += A[1];
    A[4] += A[2]*A[1];
    temp2 = 2*A[4] + LeftShift(sqr(A[5]),1);
    XGCD(g, a, z, sonsgh, temp2);
    //cout<<g<<endl;
  }

  J[0] = sqr(A[0])/(sgsh);
  J[1] = sqr(A[1])*sg/power(sg1,3);
  J[3] = sh/A[3];

  XGCD(g, a, b, sa = sh*sg1, sb = sqr(umod));
  u = J[2] = (A[2] - y*(power(A[2],3)+F))*a*sa % (sa*sb);
  
  temp = LeftShift(power(A[5],3),1) - G;
  XGCD(g, a, b, sa = sg/sg1, sb = sqr(A[1]/sg1));
  w = J[5] = (A[5] - z*temp)*a*sa % (sa*sb);

  XGCD(g, a, b, sa = sg*A[3], sonsgh2);
  J[4] = (A[4] + J[2]*(J[5]-A[5]) + z*(J[2]*temp + LeftShift(2*G*A[5],1) - A[4]*(temp2-A[4])))*a*sa % (sa*sonsgh2);

  // "Reduce."
  
  J[2] %= (J[0]/J[1]);
  J[5] %= J[1];
  J[4] += u*(J[5]-w);
  J[4] %= (J[0]/J[3]);

  //  cout<<J[0]<<J[1]<<J[2]<<J[3]<<J[4]<<J[5]<<endl;

  return(deg(sgsh)-deg(sg1)-deg(A[3]));

}

// Finds J = A*B in the case that J is primitive.

int primcomp(vec_ZZ_pX A, vec_ZZ_pX B, vec_ZZ_pX P, int g, vec_ZZ_pX &J) {

  ZZ_pX sga, sha, sgb, shb; 
  ZZ_pX d, d1, dG, dH, dp, d2; 
  ZZ_pX r1, r2, r3, r4; 
  ZZ_pX w;

  //cout<<A[0]<<A[1]<<A[2]<<A[3]<<A[4]<<A[5]<<endl;
  //cout<<B[0]<<B[1]<<B[2]<<B[3]<<B[4]<<B[5]<<endl;
  //cout<<P[0]<<P[1]<<P[2]<<P[3]<<P[4]<<P[5]<<endl;

  GCD(sga, A[0], G);
  //  GCD(sha, A[0], H);
  if(IsZero(ConstTerm(A[0])))
    SetX(sha);
  else
    sha = ZZ_pX(0,1);

  GCD(sgb, B[0], G);
  //  GCD(shb, B[0], H);
  if(IsZero(ConstTerm(B[0])))
    SetX(shb);
  else
    shb = ZZ_pX(0,1);

  if(IsOne(sha) || IsOne(shb)){
    dp = ZZ_pX(0,1);
    dH = ZZ_pX(0,1);
  }
  else{
    SetX(dp);
    if(IsOne(sha/A[3]) || IsOne(shb/B[3]))
      dH = ZZ_pX(0,1);
    else
      SetX(dH);
  }

  XGCD(d, r1, r2, A[0]/(A[1]*dp), B[0]/(B[1]*dp));

  GCD(dG, d, G); 
  GCD(d1, d, A[2]-B[2]);
  d1/=dG;

  r2 = d/d1;

  // Set s, s', and s".
  J[0] = A[0]*B[0]/(r2*dH);
  J[1] = A[1]*B[1]*r2;
  J[3] = A[3]*B[3]*dH;

  //cout<<"s,s',s'' = "<<J[0]<<J[1]<<J[3]<<endl;

  // Set u.
  J[2] = A[2] + (B[2]-A[2])*(A[0]*r1/(A[1]*d*dp));

  XGCD(d2, r3, r4, 3*sqr(J[2]), d1*sqr(dp)/dH);
  J[2] -= r3*(power(J[2], 3) + LeftShift(G, 2))/d2;
  J[2] %= (J[0]/J[1]);

  // Set v and w.

  J[4] = P[0]*A[3]*B[3]*(A[4]*B[4]+(A[5]+B[5])*LeftShift(G, 1)) + P[1]*(B[1]*A[3]*(B[2]*A[4]+LeftShift(G,1))) + P[2]*B[0]*A[3]*A[4] + P[3]*A[1]*(A[2]*B[4]+LeftShift(G,1)) + P[4]*A[1]*B[1]*A[2]*B[2] + P[5]*A[0]*B[3]*B[4];

  w = P[0]*A[3]*B[3]*(A[4]*B[5]+B[4]*A[5]+G) + P[1]*B[1]*A[3]*(B[2]*A[5]+A[4]) + P[2]*B[0]*A[3]*A[5] + P[3]*A[1]*(A[2]*B[5]+B[4]) + P[4]*A[1]*B[1]*(A[2]+B[2]) + P[5]*A[0]*B[3]*B[5];

  if(g){
    if(!IsOne(J[1])){
      SetX(r2);
      InvMod(r1, r2, J[1]); 
      //cout<<J[4]<<w<<endl;
      if(coeff(J[4],0)==0)
	J[4] = RightShift(J[4], 1);
      else
	J[4]*=r1;
      if(coeff(w,0)==0)
	w = RightShift(w, 1);
      else
	w*=r1;
    }
    else{
      J[4] = RightShift(J[4], 1);
      w = RightShift(w, 1);
    }
    //cout<<J[4]<<w<<endl;
  }

  J[5] = w%J[1];
  J[4] += J[2]*(J[5]-w);
  J[4] %= J[0]/J[3];
    
  //cout<<"The product: "<<endl;
  //cout<<J[0]<<J[1]<<J[2]<<J[3]<<J[4]<<J[5]<<endl;    

  return 0;
 }

// This function divides out most, if not all, of the nonprimitive
// part of the product A*B. The function changes A and B and
// defines I as the extra factor that puts back any extra primes
// that were divided out.
// Returns the degree of the factor of H(x) that was divided out.

int divprimes(vec_ZZ_pX &A, vec_ZZ_pX &B, ZZ_pX g, vec_ZZ_pX &I){

  ZZ_pX sha, shb, sga, sgb, sgap, sgbp;
  ZZ_pX sta, stb, stap, stbp;
  ZZ_pX d1, d1t, d2, d2t, d3, d3t, d1H, d2H, dH, dHpp, dGp, d1G, d2G;
  ZZ_pX d123, u, w;
  ZZ_pX g2, r1, r2;

  //cout<<A[0]<<A[1]<<A[2]<<A[3]<<A[4]<<A[5]<<endl;  
  //cout<<B[0]<<B[1]<<B[2]<<B[3]<<B[4]<<B[5]<<endl;  
  
  if(IsZero(ConstTerm(A[0])))
    SetX(sha);
  else
    sha = ZZ_pX(0,1);
  if(IsZero(ConstTerm(B[0])))
    SetX(shb);
  else
    shb = ZZ_pX(0,1);

  GCD(sga, A[0], G);
  GCD(sgb, B[0], G);
  GCD(sgap, A[1], G);
  GCD(sgbp, B[1], G);
  
  GCD(dGp, sgap, sgbp);
  if(IsOne(A[3]) || IsOne(B[3]))
    dHpp = ZZ_pX(0, 1);
  else
    SetX(dHpp);

  //GCD(dHpp, A[3], B[3]);
  if(!IsOne(sha) && !IsOne(shb)) 
    dH = GCD(sha, shb)/GCD(sha/A[3], shb/B[3]);
  else
    dH = ZZ_pX(0, 1);

  sta = A[0]/(sga*sha);
  stb = B[0]/(sgb*shb);
  stap = A[1]/sgap;
  stbp = B[1]/sgbp;

  GCD(d1H, sha/A[3], B[3]);
  GCD(d2H, shb/B[3], A[3]);
  GCD(d1G, sga/sgap, sgbp);
  GCD(d2G, sgb/sgbp, sgap);
  
  d1t =GCD(GCD(sta/stap, stbp), A[2]+LeftShift(B[5],1));
  d1 = d1G*d1t;
  d2t = GCD(GCD(stb/stbp, stap), B[2]+LeftShift(A[5], 1));
  d2 = d2G*d2t;
  d3t = GCD(GCD(sta/d2t, stb/d1t), A[4]+B[4]+LeftShift(A[5]*B[5], 1));
  d3 = dGp*d3t;
  d123=d1*d2*d3;

  // Create I;

  I[0] = dHpp*d3;
  I[1] = ZZ_pX(0,1);

  if(IsOne(d3t))
    I[2] = ZZ_pX();
  else{  
    I[2] = LeftShift(A[5]+B[5], 1);
    if(!IsOne(dHpp))
      I[2]*=LeftShift(InvMod(dHpp, d3t),1);
    I[2] %= I[0];
  }
  I[3] = ZZ_pX(0,1);
  I[4] = LeftShift(-A[5]*B[5], 1)%I[0];
  I[5] = ZZ_pX();
    
  A[0]/=(d123*dH);
  A[1]/=(d2*d3);
  A[3]/=(d2H*dHpp);

  B[0]/=(d123*dH);
  B[1]/=(d1*d3);
  B[3]/=(d1H*dHpp);

  // "Reduce" the coefficients u, v, and w in A and B
  u = A[2];
  w = A[5];
  A[2] %= (A[0]/A[1]);
  A[5] %= A[1];
  A[4] += u*(A[5]-w);
  A[4] %= (A[0]/A[3]);

  u = B[2];
  w = B[5];
  B[2] %= (B[0]/B[1]);
  B[5] %= B[1];
  B[4] += u*(B[5]-w);
  B[4] %= (B[0]/B[3]);

  //cout<<"New A, B, and I:"<<endl;
  //cout<<A[0]<<A[1]<<A[2]<<A[3]<<A[4]<<A[5]<<endl;  
  //cout<<B[0]<<B[1]<<B[2]<<B[3]<<B[4]<<B[5]<<endl;  
  //cout<<I[0]<<I[1]<<I[2]<<I[3]<<I[4]<<I[5]<<endl;

  return deg(dH);
}

void babystep0(ZZ_p half, ZZ_p U){
  ZZ_pX a0, a1, a2; // neighbor
  ZZ_pX o0, o1, o2, i0, i1, i2; // temp variables
  int i;
  ZZ_pX dold;  //old denominator
  ZZ_pX autnu;
  ZZ_pX g;            // gcd of numerator and denominator

  // Step 2 is reducing {1, rho, omega}.
  reduce(half);
  // Step 3.
  // Step 3.1
  
  // Step 3.2
  //  if(deg(autnu = aut(nu0, nu1, nu2, U)) == ddeg){
  if(deg(autnu = LeftShift(nu0, shift) + nu1*U*rho + nu2*U*U*omega) == (shift+ddeg)){
    
    a0 = nu0-d*ZZ_p(LeadCoeff(autnu));
    a1 = nu1; a2=nu2;

    // Update the degree.
    d0 += eltdeg(a0,a1,a2);
    d1 += eltdeg(a0, U*a1, U*U*a2);
    d2 += eltdeg(a0, U*U*a1, U*a2);

    //cout<<"1. "<<eltdeg(a0,a1,a2)<<" "<<eltdeg(nu0, nu1, nu2)<<endl;
    cout<<"Stepping in the direction: ("<<eltdeg(a0,a1,a2)<<", "<<eltdeg(a0, U*a1, U*U*a2)<<", "<<eltdeg(a0, U*U*a1, U*a2)<<")."<<endl;

    //Step 3.3
    // theta := theta*alpha
    div(o0, LeftShift((t1*a2 + t2*a1)*G,1) + t0*a0, d);
    div(o1, t0*a1 + t1*a0 + t2*a2*G, d);
    div(o2, t0*a2 + LeftShift(t1*a1,1) + t2*a0, d);
    t0 = o0; t1 = o1; t2 = o2;
     
    // Step 3.4
    // mu = a^-1 nu = mu*a^-1
    // a^-1
    o0 = (sqr(a0) + LeftShift(-a1*a2*G,1));
    o1 = (sqr(a2)*G - a0*a1);
    o2 = (-a0*a2 + LeftShift(sqr(a1),1));
    
    i0 = LeftShift((o1*mu2 + o2*mu1)*G,1) + o0*mu0;
    i1 = o0*mu1 + o1*mu0 + o2*mu2*G;
    i2 = o0*mu2 + LeftShift(o1*mu1,1) + o2*mu0;
    
    nu0 = i0; nu1 = i1; nu2 = i2;
    mu0 = o0*d; mu1 = o1*d; mu2 = o2*d;
    
    // Set the denominator
    d=LeftShift(a1*sqr(a1)*G,2);
    d+=LeftShift(sqr(a2*G)*a2-3*a0*a1*a2*G,1);
    d+=a0*sqr(a0);
     
  }
  else{
    // Step 3.3
    // Update the degree.
    d0+=eltdeg(mu0, mu1,  mu2);
    d1+=eltdeg(mu0, U*mu1, U*U*mu2);
    d2+=eltdeg(mu0, U*U*mu1, U*mu2);

    //cout<<"2. "<<eltdeg(mu0,mu1,mu2)<<endl;
    cout<<"Stepping in the direction: ("<<eltdeg(mu0,mu1,mu2)<<", "<<eltdeg(mu0, U*mu1, U*U*mu2)<<", "<<eltdeg(mu0, U*U*mu1, U*mu2)<<")."<<endl;
    
    // theta := theta*mu
    div(o0, LeftShift((t1*mu2 + t2*mu1)*G,1) + t0*mu0, d);
    div(o1, t0*mu1 + t1*mu0 + t2*mu2*G, d);
    div(o2, t0*mu2 + LeftShift(t1*mu1,1) + t2*mu0, d);
    t0 = o0; t1 = o1; t2 = o2;
     
    // Step 3.4
    // mu = mu^-1 nu = nu*mu^-1
    // Set the denominator.
    dold=d;
    d=LeftShift(mu1*sqr(mu1)*G,2);
    d+=LeftShift(sqr(mu2*G)*mu2-3*mu0*mu1*mu2*G,1);
    d+=mu0*sqr(mu0);
    
    // mu^-1
    o0 = (sqr(mu0) + LeftShift(-mu1*mu2*G,1));
    o1 = (sqr(mu2)*G - mu0*mu1);
    o2 = (-mu0*mu2 + LeftShift(sqr(mu1),1));
    
    i0 = LeftShift((o1*nu2 + o2*nu1)*G,1) + o0*nu0;
    i1 = o0*nu1 + o1*nu0 + o2*nu2*G;
    i2 = o0*nu2 + LeftShift(o1*nu1,1) + o2*nu0;
    
    nu0 = i0; nu1 = i1; nu2 = i2;
    mu0 = o0*dold; mu1 = o1*dold; mu2 = o2*dold;  
  }
  
  // Check to make sure gcd(d, mu, nu) == 1.
  g = GCD(d, nu2);
  if(g!=1){
    g = GCD(g, nu1); g = GCD(g, nu0);
    if(g!=1){
      g = GCD(g, mu2); g = GCD(g, mu1); g = GCD(g, mu0);
      if(g!=1) {
        d/=g; mu0/=g; mu1/=g; mu2/=g; nu0/=g; nu1/=g; nu2/=g;
      }
    }
  }
  ddeg=deg(d);
  
  //    cout<<"gcd = "<<g<<endl;
  
  // Step 3.5 Reduce the ideal {1, mu, nu}.
  reduce(half);
  
  //Normalize d;
  d*=inv(LeadCoeff(d));
  // cout<<"MU_{("<<l<<", "<<m<<")} = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
  // cout<<"NU_{("<<l<<", "<<m<<")} = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";


  return;
  
}

// Returns the degree of the step.

void babystep1(ZZ_p half, ZZ_p U){
  int cont = 1;
  int j;
  ZZ_pX dold;
  ZZ_pX o0, o1, o2, i0, i1, i2; // temp variables
  ZZ_pX b0, b1, b2; // beta
  ZZ_pX tau;
  ZZ_pX g; // gcd
  //p=0;

  // Step 2.1
  reduce1(half, U);
  
  // Step 2.2
  if(deg(tau = LeftShift(tau0,shift)+tau1*U*U*rho+tau2*U*omega) == ddeg+shift){
    b0 = tau0 - LeadCoeff( tau )*d; b1 = tau1; b2 = tau2;

    // Update the degrees.
    d0 += eltdeg(b0, b1, b2);
    d1 += eltdeg(b0, U*b1, U*U*b2);
    d2 += eltdeg(b0, U*U*b1, U*b2);

    cout<<"Stepping in the direction: ("<<eltdeg(b0,b1,b2)<<", "<<eltdeg(b0, U*b1, U*U*b2)<<", "<<eltdeg(b0, U*U*b1, U*b2)<<")."<<endl;
    
    //Step 2.3
    // phi := phi*beta
    div(o0, LeftShift((phi1*b2 + phi2*b1)*G,1) + phi0*b0, d);
    div(o1, phi0*b1 + phi1*b0 + phi2*b2*G, d);
    div(o2, phi0*b2 + LeftShift(phi1*b1,1) + phi2*b0, d);
    phi0 = o0; phi1 = o1; phi2 = o2;
    
    //Step 2.4
    // sigma = beta^-1 tau = sigma*beta^-1
    
    // beta^-1
    o0 = (sqr(b0) + LeftShift(-b1*b2*G,1));
    o1 = (sqr(b2)*G - b0*b1);
    o2 = (-b0*b2 + LeftShift(sqr(b1),1));
    
    i0 = LeftShift((o1*sig2 + o2*sig1)*G,1) + o0*sig0;
    i1 = o0*sig1 + o1*sig0 + o2*sig2*G;
    i2 = o0*sig2 + LeftShift(o1*sig1,1) + o2*sig0;
    
    tau0 = i0; tau1 = i1; tau2 = i2;
    sig0 = o0*d; sig1 = o1*d; sig2 = o2*d;
    
    // Set the denominator
    d=LeftShift(b1*sqr(b1)*G,2);
    d+=LeftShift(sqr(b2*G)*b2-3*b0*b1*b2*G,1);
    d+=b0*sqr(b0);
    //Normalize d;
    d*=inv(LeadCoeff(d));
    
    // Check to make sure gcd(d, sig, tau) == 1.
    g = GCD(d, tau2);
    if(g!=1){
      g = GCD(g, tau1); g = GCD(g, tau0);
      if(g!=1){
        g = GCD(g, sig2); g = GCD(g, sig1); g = GCD(g, sig0);
        if(g!=1) {
          d/=g; sig0/=g; sig1/=g; sig2/=g; tau0/=g; tau1/=g; tau2/=g;
        }
      }
    }
    
    ddeg=deg(d);      
    
  }
  else{
    //Step 2.3

    // Update the degrees.
    d0 += eltdeg(sig0, sig1, sig2);
    d1 += eltdeg(sig0, U*sig1, U*U*sig2);
    d2 += eltdeg(sig0, U*U*sig1, U*sig2);
    
    cout<<"Stepping in the direction: ("<<eltdeg(sig0,sig1,sig2)<<", "<<eltdeg(sig0, U*sig1, U*U*sig2)<<", "<<eltdeg(sig0, U*U*sig1, U*sig2)<<")."<<endl;

    // phi := phi*sigma
    div(o0, LeftShift((phi1*sig2 + phi2*sig1)*G,1) + phi0*sig0, d);
    div(o1, phi0*sig1 + phi1*sig0 + phi2*sig2*G, d);
    div(o2, phi0*sig2 + LeftShift(phi1*sig1,1) + phi2*sig0, d);
    phi0 = o0; phi1 = o1; phi2 = o2;
    
    //Step 2.4
    // sigma = sigma^-1 tau = tau*sigma^-1
    
    // Set the denominator
    dold=d;
    d=LeftShift(sig1*sqr(sig1)*G,2);
    d+=LeftShift(sqr(sig2*G)*sig2-3*sig0*sig1*sig2*G,1);
    d+=sig0*sqr(sig0);
    
    // sigma^-1
    o0 = (sqr(sig0) + LeftShift(-sig1*sig2*G,1));
    o1 = (sqr(sig2)*G - sig0*sig1);
    o2 = (-sig0*sig2 + LeftShift(sqr(sig1),1));
    
    i0 = LeftShift((o1*tau2 + o2*tau1)*G,1) + o0*tau0;
    i1 = o0*tau1 + o1*tau0 + o2*tau2*G;
    i2 = o0*tau2 + LeftShift(o1*tau1,1) + o2*tau0;
    
    tau0 = i0; tau1 = i1; tau2 = i2;
    sig0 = o0*dold; sig1 = o1*dold; sig2 = o2*dold;
    
    
    //Normalize d;
    d*=inv(LeadCoeff(d));
    
    // Check to make sure gcd(d, sig, tau) == 1.
    g = GCD(d, tau2);
    if(g!=1){
      g = GCD(g, tau1); g = GCD(g, tau0);
      if(g!=1){
        g = GCD(g, sig2); g = GCD(g, sig1); g = GCD(g, sig0);
        if(g!=1) {
          d/=g; sig0/=g; sig1/=g; sig2/=g; tau0/=g; tau1/=g; tau2/=g;
        }
      }
    }
    
    ddeg=deg(d);     
  }
  
  //Step 2.5
  mu0=sig0; mu1=sig1; mu2=sig2;
  nu0=tau0; nu1=tau1; nu2=tau2;
  
  reduce(half);

  return;
}

// Returns the degree of the step.

void babystep2(ZZ_p half, ZZ_p U){
  int cont = 1;
  int j;
  ZZ_pX dold;
  ZZ_pX o0, o1, o2, i0, i1, i2; // temp variables
  ZZ_pX b0, b1, b2; // beta
  ZZ_pX tau;
  ZZ_pX g; // gcd
  //p=0;

  // Step 2.1
  reduce2(half, U);
  
  // Step 2.2
  if(deg(tau = LeftShift(tau0,shift)+tau1*rho+tau2*omega) == ddeg+shift){
    b0 = tau0 - LeadCoeff( tau )*d; b1 = tau1; b2 = tau2;

    //cout<<"("<<d0<<", "<<d1<<", "<<d2<<")"<<endl;

    // Update the degrees.
    d0 += eltdeg(b0, b1, b2);
    d1 += eltdeg(b0, U*b1, U*U*b2);
    d2 += eltdeg(b0, U*U*b1, U*b2);
    //cout<<"("<<d0<<", "<<d1<<", "<<d2<<")"<<endl;
    cout<<"Stepping in the direction: ("<<eltdeg(b0,b1,b2)<<", "<<eltdeg(b0, U*b1, U*U*b2)<<", "<<eltdeg(b0, U*U*b1, U*b2)<<")."<<endl;
    
    //Step 2.3
    // phi := phi*beta
    div(o0, LeftShift((phi1*b2 + phi2*b1)*G,1) + phi0*b0, d);
    div(o1, phi0*b1 + phi1*b0 + phi2*b2*G, d);
    div(o2, phi0*b2 + LeftShift(phi1*b1,1) + phi2*b0, d);
    phi0 = o0; phi1 = o1; phi2 = o2;
    
    //Step 2.4
    // sigma = beta^-1 tau = sigma*beta^-1
    
    // beta^-1
    o0 = (sqr(b0) + LeftShift(-b1*b2*G,1));
    o1 = (sqr(b2)*G - b0*b1);
    o2 = (-b0*b2 + LeftShift(sqr(b1),1));
    
    i0 = LeftShift((o1*sig2 + o2*sig1)*G,1) + o0*sig0;
    i1 = o0*sig1 + o1*sig0 + o2*sig2*G;
    i2 = o0*sig2 + LeftShift(o1*sig1,1) + o2*sig0;
    
    tau0 = i0; tau1 = i1; tau2 = i2;
    sig0 = o0*d; sig1 = o1*d; sig2 = o2*d;
    
    // Set the denominator
    d=LeftShift(b1*sqr(b1)*G,2);
    d+=LeftShift(sqr(b2*G)*b2-3*b0*b1*b2*G,1);
    d+=b0*sqr(b0);
    //Normalize d;
    d*=inv(LeadCoeff(d));
    
    // Check to make sure gcd(d, sig, tau) == 1.
    g = GCD(d, tau2);
    if(g!=1){
      g = GCD(g, tau1); g = GCD(g, tau0);
      if(g!=1){
        g = GCD(g, sig2); g = GCD(g, sig1); g = GCD(g, sig0);
        if(g!=1) {
          d/=g; sig0/=g; sig1/=g; sig2/=g; tau0/=g; tau1/=g; tau2/=g;
        }
      }
    }
    
    ddeg=deg(d);      
    
  }
  else{
    //Step 2.3
//cout<<"("<<d0<<", "<<d1<<", "<<d2<<")"<<endl;
    // Update the degrees.
    d0 += eltdeg(sig0, sig1, sig2);
    d1 += eltdeg(sig0, U*sig1, U*U*sig2);
    d2 += eltdeg(sig0, U*U*sig1, U*sig2);
//cout<<"("<<d0<<", "<<d1<<", "<<d2<<")"<<endl;
    cout<<"Stepping in the direction: ("<<eltdeg(sig0,sig1,sig2)<<", "<<eltdeg(sig0, U*sig1, U*U*sig2)<<", "<<eltdeg(sig0, U*U*sig1, U*sig2)<<")."<<endl;

    // phi := phi*sigma
    div(o0, LeftShift((phi1*sig2 + phi2*sig1)*G,1) + phi0*sig0, d);
    div(o1, phi0*sig1 + phi1*sig0 + phi2*sig2*G, d);
    div(o2, phi0*sig2 + LeftShift(phi1*sig1,1) + phi2*sig0, d);
    phi0 = o0; phi1 = o1; phi2 = o2;
    
    //Step 2.4
    // sigma = sigma^-1 tau = tau*sigma^-1
    
    // Set the denominator
    dold=d;
    d=LeftShift(sig1*sqr(sig1)*G,2);
    d+=LeftShift(sqr(sig2*G)*sig2-3*sig0*sig1*sig2*G,1);
    d+=sig0*sqr(sig0);
    
    // sigma^-1
    o0 = (sqr(sig0) + LeftShift(-sig1*sig2*G,1));
    o1 = (sqr(sig2)*G - sig0*sig1);
    o2 = (-sig0*sig2 + LeftShift(sqr(sig1),1));
    
    i0 = LeftShift((o1*tau2 + o2*tau1)*G,1) + o0*tau0;
    i1 = o0*tau1 + o1*tau0 + o2*tau2*G;
    i2 = o0*tau2 + LeftShift(o1*tau1,1) + o2*tau0;
    
    tau0 = i0; tau1 = i1; tau2 = i2;
    sig0 = o0*dold; sig1 = o1*dold; sig2 = o2*dold;
    
    
    //Normalize d;
    d*=inv(LeadCoeff(d));
    
    // Check to make sure gcd(d, sig, tau) == 1.
    g = GCD(d, tau2);
    if(g!=1){
      g = GCD(g, tau1); g = GCD(g, tau0);
      if(g!=1){
        g = GCD(g, sig2); g = GCD(g, sig1); g = GCD(g, sig0);
        if(g!=1) {
          d/=g; sig0/=g; sig1/=g; sig2/=g; tau0/=g; tau1/=g; tau2/=g;
        }
      }
    }
    
    ddeg=deg(d);     
  }
  
  //Step 2.5
  mu0=sig0; mu1=sig1; mu2=sig2;
  nu0=tau0; nu1=tau1; nu2=tau2;
  
  reduce(half);

  return;
}

// Searches all ideals to see if one of them is the current (1, mu, nu).
// Searches the range lo to hi.
// Returns -1 if there is no match.
// Returns the location if there is.

int search(int lo, int hi, int verbose){

  int i;
  
  for(i=lo; i< hi; i++){
    if(D[i]==d){
      // Now check for mu and nu.
      if((I[6*i] == mu0) && (I[6*i+1] == mu1) && (I[6*i+2] == mu2) && (I[6*i+3] == nu0) && (I[6*i+4] == nu1) && (I[6*i+5] == nu2)){
	if(verbose){
	  cout<<"MATCH = ("<<coordl[i]<<", "<<coordm[i]<<") = \n(";
	  cout<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
	  cout<<"("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n";	
	  cout<<"Degree = ("<<deg0[i]<<", "<<deg1[i]<<", "<<deg2[i]<<").\n\n";
	}
	return(i);
      }
    }
  }

  // If it gets here, there's no ideal in the list.
  return (-1);

}

// Inserts an ideal into the list, identifies it with location (l0, m0).

void insert(int l0, int m0, int verbose){

  if(verbose){
    cout<<"MU_{("<<l0<<", "<<m0<<")} = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
    cout<<"NU_{("<<l0<<", "<<m0<<")} = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n";
    cout<<"Degree = ("<<d0<<", "<<d1<<", "<<d2<<").\n\n";
  }
  //C[3*num] = phi0; C[3*num+1]=phi1; C[3*num+2]=phi2;
  I[6*num] = mu0; I[6*num+1]=mu1; I[6*num+2]=mu2;
  I[6*num+3] = nu0; I[6*num+4]=nu1; I[6*num+5]=nu2;
  D[num] = d;
  coordl[num] = l0;
  coordm[num] = m0;
  deg0[num] = d0;
  deg1[num] = d1;
  deg2[num] = d2;
  num++;


  // ddeg=deg(d);

}

// Input a reduced basis of an ideal.
// Output a reduced basis of an equivalent reduced ideal.
// Return the degree \deg(\psi), where $f_{new} = (\psi^{-1})f_{old}
// Algorithm 6.11 in Ideal Arithmetic and Infrastructure 
// in Purely Cubic Function Fields, Scheidler.

void reduceideal(ZZ_p half, ZZ_p U){
  ZZ_pX dold, o0, o1, o2, i0, i1, i2, g;
  ZZ_pX a0,a1,a2;
  //int psideg=0; // The degree.
  int degtemp, ed;

  cout<<"("<<d0<<", "<<d1<<", "<<d2<<")"<<endl;

  // Step 2
  while( (degtemp = eltdeg(mu0,mu1,mu2) ) <= 0 ){
    //psideg += degtemp;

    //cout<<"First loop: "<<xideg(nu1, nu2)<<" < "<<xideg(mu1, mu2)<<"; "<<etadeg(mu1, mu2)<<" < 0 <= "<<etadeg(nu1, nu2)<<"; "<<zetadeg(mu0, mu1, mu2)<<" < 0; "<<zetadeg(nu0, nu1, nu2)<<" < 0."<<endl;

    // Update the degree.
    d0 += degtemp;
    d1 += eltdeg(mu0, U*mu1, U*U*mu2);
    d2 += eltdeg(mu0, U*U*mu1, U*mu2);
    
    cout<<"N(mu) = "<<sqr(mu0)*mu0 + LeftShift(sqr(mu1)*mu1*G, 2) + LeftShift(sqr(mu2)*mu2*sqr(G), 1) - LeftShift(3*mu0*mu1*mu2*G, 1)<<endl;
    cout<<"Step 1: ("<<d0<<", "<<d1<<", "<<d2<<")"<<endl;
    // Step 2.1
    // Set (mu, nu) = (mu^{-1}, nu*mu^{-1})
    dold=d;
    d=LeftShift(mu1*sqr(mu1)*G,2);
    d+=LeftShift(sqr(mu2*G)*mu2-3*mu0*mu1*mu2*G,1);
    d+=mu0*sqr(mu0);  

    // mu^-1
    o0 = (sqr(mu0) + LeftShift(-mu1*mu2*G,1));
    o1 = (sqr(mu2)*G - mu0*mu1);
    o2 = (-mu0*mu2 + LeftShift(sqr(mu1),1));
    
    i0 = LeftShift((o1*nu2 + o2*nu1)*G,1) + o0*nu0;
    i1 = o0*nu1 + o1*nu0 + o2*nu2*G;
    i2 = o0*nu2 + LeftShift(o1*nu1,1) + o2*nu0;
    
    nu0 = i0; nu1 = i1; nu2 = i2;
    mu0 = o0*dold; mu1 = o1*dold; mu2 = o2*dold;

    // Check to make sure gcd(d, mu, nu) == 1.
    g = GCD(d, nu2);
    if(g!=1){
      g = GCD(g, nu1); g = GCD(g, nu0);
      if(g!=1){
	g = GCD(g, mu2); g = GCD(g, mu1); g = GCD(g, mu0);
	if(g!=1) {
	  d/=g; mu0/=g; mu1/=g; mu2/=g; nu0/=g; nu1/=g; nu2/=g;
	}
      }
    }
    ddeg = deg(d);
    
    // Normalize d.
    d*=inv(LeadCoeff(d));    

    // Step 2.2
    reduce(half);    
  }

  // Somehow fix nu if |nu| = |eta_nu| = 0.
  //ed=etadeg(nu1, nu2);
  //g = LeftShift(nu0, shift) + nu1*rho + nu2*omega;
  //degtemp = deg(g)-shift-ddeg;

  // cout<<"etadeg = "<<ed<<", |nu| = "<<degtemp<<endl;

  //if((degtemp == ed) && (ed == 0)){
  //  nu0 -= d*LeadCoeff(g);
  // }

  //ed=etadeg(nu1, nu2);
  //g = LeftShift(nu0, shift) + nu1*rho + nu2*omega;
  //degtemp = deg(g)-shift-ddeg;

  // cout<<"etadeg = "<<ed<<", |nu| = "<<degtemp<<endl;
  // Step 3.
  // If |nu| < |eta_{nu}| = 1, then set (mu, nu) = (mu*nu^{-1}, nu^{-1})
  // If |nu| < |eta_{nu}| = 1, then set (mu, nu) = (mu*a^{-1}, a^{-1}),
  // where a = alpha = nu - sgn(nu^{(1)}).
  

  if( ((degtemp = eltdeg(nu0, nu1, nu2)) <= 0) && ((ed=etadeg(nu1, nu2)) <= 0) ){
    //cout<<xideg(nu1, nu2)<<" < "<<xideg(mu1, mu2)<<"; "<<etadeg(mu1, mu2)<<" < 0 <= "<<etadeg(nu1, nu2)<<"; "<<zetadeg(mu0, mu1, mu2)<<" < 0; "<<zetadeg(nu0, nu1, nu2)<<" < 0."<<endl;
    
    a0 = nu0 - d*ZZ_p(LeadCoeff(aut(nu0, nu1, nu2, U)));
    a1 = nu1; a2=nu2;

    //psideg += degtemp;
    // Update the degree.
    d0 += eltdeg(a0,a1,a2);
    d1 += eltdeg(a0, U*a1, U*U*a2);
    d2 += eltdeg(a0, U*U*a1, U*a2);
    cout<<"N(alpha) = "<<sqr(a0)*a0 + LeftShift(sqr(a1)*a1*G, 2) + LeftShift(sqr(a2)*a2*sqr(G), 1) - LeftShift(3*a0*a1*a2*G, 1)<<endl;
    cout<<"Step 2: ("<<d0<<", "<<d1<<", "<<d2<<")"<<endl;
    dold=d;
    d=LeftShift(a1*sqr(a1)*G,2);
    d+=LeftShift(sqr(a2*G)*a2-3*a0*a1*a2*G,1);
    d+=a0*sqr(a0);  

    // alpha^-1
    o0 = (sqr(a0) + LeftShift(-a1*a2*G,1));
    o1 = (sqr(a2)*G - a0*a1);
    o2 = (-a0*a2 + LeftShift(sqr(a1),1));
    
    i0 = LeftShift((o1*mu2 + o2*mu1)*G,1) + o0*mu0;
    i1 = o0*mu1 + o1*mu0 + o2*mu2*G;
    i2 = o0*mu2 + LeftShift(o1*mu1,1) + o2*mu0;
     
    // Set (mu, nu) = (mu*alpha^{-1}, alpha^{-1})
    mu0 = i0; mu1 = i1; mu2 = i2;
    nu0 = o0*dold; nu1 = o1*dold; nu2 = o2*dold;

    // Check to make sure gcd(d, mu, nu) == 1.
    g = GCD(d, nu2);
    if(g!=1){
      g = GCD(g, nu1); g = GCD(g, nu0);
      if(g!=1){
	g = GCD(g, mu2); g = GCD(g, mu1); g = GCD(g, mu0);
	if(g!=1) {
	  d/=g; mu0/=g; mu1/=g; mu2/=g; nu0/=g; nu1/=g; nu2/=g;
	}
      }
    }
    ddeg = deg(d);
    
    // Normalize d.
    d*=inv(LeadCoeff(d));  

    //cout<<xideg(nu1, nu2)<<" < "<<xideg(mu1, mu2)<<"; "<<etadeg(mu1, mu2)<<" < 0 <= "<<etadeg(nu1, nu2)<<"; "<<zetadeg(mu0, mu1, mu2)<<" < 0; "<<zetadeg(nu0, nu1, nu2)<<" < 0."<<endl;

    // Step 3.2
    reduce(half);

    //cout<<xideg(nu1, nu2)<<" < "<<xideg(mu1, mu2)<<"; "<<etadeg(mu1, mu2)<<" < 0 <= "<<etadeg(nu1, nu2)<<"; "<<zetadeg(mu0, mu1, mu2)<<" < 0; "<<zetadeg(nu0, nu1, nu2)<<" < 0."<<endl;   

    ed=etadeg(nu1, nu2);
    g = LeftShift(nu0, shift) + nu1*rho + nu2*omega;
    degtemp = deg(g)-shift-ddeg;
  
    //cout<<"deg(mu) = "<<eltdeg(mu0, mu1, mu2)<<" etadeg = "<<ed<<", deg(nu) = "<<degtemp<<endl;

  }

  return;
}

int main(){

  int i, j;
  ZZ_pX autnu;
  int coefg[5]; //= {1,1,0,2,1};
  int badinput = 1;
  long ilen, dlen; // clen;
  int cont = 1;
  int step;
  int match;
  vec_ZZ_pX A, B, J;
  int c1, c2;
  int divtemp;

  J.SetLength(6);
  A.SetLength(7);
  B.SetLength(6);

  while(badinput){
    cout<<"Characteristic of the base field, P: ";
    cin>>P;

    if(P%3 == 1)
      badinput = 0;
    else
      cout<<"The characteristic must be P = 1 (mod 3).\n";
  }

  // This needs to change.
  dlen = P*P*P + 6*P*P*SqrRoot(P) + 15*P*P + 20*P*SqrRoot(P) + 15*P + 6*SqrRoot(P) + 1;
  ilen = 6*dlen;
  //clen = 3*dlen;

  //C.SetLength(clen);
  I.SetLength(ilen);
  D.SetLength(dlen);
  
  cout<<"Coefficients of G = g_0 + g_1*x + g_2*x^2 + g_3*x^3 + x^4:\n";
  for( i=0; i<4; i++){
    cout<<"g_"<<i<<": ";
    cin>>coefg[i];
  }

  // Set the characterisitc.
  ZZ_p::init(to_ZZ(P));

  ZZ_p half = inv(to_ZZ_p(2));


  // Find a primitive cube root of 1
  for(i=1; i<P; i++){
    if( (i*i + i + 1)%P == 0 ) {
      u=i;
      break;
    }
  }
  ZZ_p U = to_ZZ_p(u);

  // The coefficients of the basis elements.
  t0 = to_ZZ_pX(1);
  mu0 = ZZ_pX(0,0), mu1 = ZZ_pX(0,1), mu2 = ZZ_pX(0,0); 
  nu0 = ZZ_pX(0,0), nu1 = ZZ_pX(0,0), nu2 = ZZ_pX(0,1);

  // Set the denominator
  d = to_ZZ_pX(1);

  // Initialize rho, omega, and G.
  
  for(i=0; i<4;i++)
    G += ZZ_pX(i,coefg[i]);
  G+=ZZ_pX(4, 1);

  F = LeftShift(G,2);

  // Find a good approximation for the basis elements.
  basis(shift);

  // Run the reduction.
  reduce(half);

  //  cout<<"MU_{("<<l<<", "<<m<<")} = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
  //  cout<<"NU_{("<<l<<", "<<m<<")} = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
 
  l=0;
  cont=1;
  num = 0;
  //Get the initial period.
  do{
    insert(l, 0, 1);

    ddeg=deg(d);
    l++;

    babystep0(half, U);

    // See if {1, mu, nu} is in I.
    cont = 1+ search(0, l, 0);
    
  } while(!cont);

  dege1 = d0 - deg0[cont-1];

  // For keeping track of where the ideals off the main loop are located.
  int idealloc[l];

  // Check for the preperiod.
  if(p)
    l-=p;

  cout<<"PREPERIOD p = "<<p<<endl;
  cout<<"PERIOD #1 l = "<<l<<"\n"; 

  // Get the secondary periods
  for(i=p; i<p+l; i++){
    m=0;
    d0 = deg0[i];
    d1 = deg1[i];
    d2 = deg2[i];
    cout<<"Kick off from ("<<d0<<", "<<d1<<", "<<d2<<")"<<endl;
    //phi0 = C[3*i]; phi1 = C[3*i+1]; phi2 = C[3*i+2];
    sig0 = I[6*i]; sig1=I[6*i+1]; sig2=I[6*i+2];
    tau0=I[6*i+3]; tau1=I[6*i+4]; tau2 = I[6*i+5];
    d = D[i];
    ddeg=deg(d);

    cont=0;

    while(!cont){
      m++;
      babystep2(half, U);

      cont = 1 + search(0, num, 1);

      if(cont){
	if(i==p){
	  dege2 = d2 - deg0[cont-1] + deg0[i];
	}
	idealloc[i] = num-m+1;
	//cout<<i<<" "<<idealloc[i]<<endl;
	cout<<"PERIOD #2 m_"<<i<<" = "<<m<<" AT "<<endl; 
      }
      else{
	//cout<<num<<endl;
	insert(i, m, 1);	
        ddeg=deg(d);
      }
    }
  }

  // Now that we have mapped a number of points automatically,
  // We now go interactive to navigate through the infrastructure.
  
  cout<<"TOTAL SO FAR = "<<num<<endl;
  cout<<"DEGREE(e_1) = "<<dege1<<"\n";
  cout<<"DEGREE(e_2) = "<<dege2<<"\n\n";

  //t0 = C[0]; t1 = C[1]; t2 = C[2];
  mu0 = I[0]; mu1=I[1]; mu2=I[2];
  nu0=I[3]; nu1=I[4]; nu2 = I[5];
  d = D[0];
  ddeg=deg(d);
  d0=d1=d2=0;

  cont = 1;
  l--;
  while (cont){
    match=0;
    cout<<"3 (Quit); 4 (Restart); 5 (Comp); 6 (Square); 0, 1, 2 (Step): ";
    cin>>step;
    
    if(step == 3){
      cont = 0;
      break;
    }
    else if( step == 4){
      //t0 = C[0]; t1 = C[1]; t2 = C[2];
      mu0 = I[0]; mu1=I[1]; mu2=I[2];
      nu0=I[3]; nu1=I[4]; nu2 = I[5];
      d = D[0];
      ddeg=deg(d);
      d0=d1=d2=0;
      cout<<"At Degree = ("<<d0<<", "<<d1<<", "<<d2<<").\n";
    }
    else if (step == 0){
      babystep0(half, U);
      cout<<"At Degree = ("<<d0<<", "<<d1<<", "<<d2<<").\n";
      /*
      if(d0 >=dege1)
	d0%=dege1;
      if(d2 >=dege2)
	d2%=dege2;
      cout<<"At Degree = ("<<d0<<", "<<d2<<").\n";
      */
      if( search(0, num, 1) == -1 ){
	insert(-l, m, 1);
	m++;
      }

    }
    else if (step == 1){
      phi0=t0; phi1=t1; phi2=t2;
      sig0=mu0; sig1=mu1; sig2=mu2;
      tau0=nu0; tau1=nu1; tau2=nu2;

      babystep1(half, U);
      cout<<"At Degree = ("<<d0<<", "<<d1<<", "<<d2<<").\n";
      /*
      if(d0 >=dege1)
	d0%=dege1;
      if(d2 >=dege2)
	d2%=dege2;
      cout<<"At Degree = ("<<d0<<", "<<d2<<").\n";
      */
      if( search(0, num, 1) == -1 ){
	insert(-l, m, 1);
	m++;
      }

    }
    else if (step == 2){
      phi0=t0; phi1=t1; phi2=t2;
      sig0=mu0; sig1=mu1; sig2=mu2;
      tau0=nu0; tau1=nu1; tau2=nu2;

      babystep2(half, U);
      cout<<"At Degree = ("<<d0<<", "<<d1<<", "<<d2<<").\n";
      /*
      if(d0 >=dege1)
	d0%=dege1;
      if(d2 >=dege2)
	d2%=dege2;
      cout<<"At Degree = ("<<d0<<", "<<d2<<").\n";
      */

      if( search(0, num, 1) == -1 ){
	insert(-l, m, 1);
	m++;
      }
    }
    else if (step == 5){
      cout<<"Multiply by ideal. 0-chain coordinate: ";
      cin>>c1;
      cout<<"Multiply by ideal. 2-chain coordinate: ";
      cin>>c2;
      A[0] = d;
      A[1] = mu0;
      A[2] = mu1;
      A[3] = mu2;
      A[4] = nu0;
      A[5] = nu1;
      A[6] = nu2;

      triangle(A, B);

      cout<<"A: S = "<<B[0]<<", S' = "<<B[1]<<", S'' = "<<B[3]<<endl;

      //cout<<"MU = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
      //cout<<"NU = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";

      if(c2>0){
	i = idealloc[c1]+c2-1;
	//cout<<idealloc[c1]<<endl;
      }
      else
	i = c1;

      A[0] = D[i];
      A[1] = I[6*i];
      A[2] = I[6*i+1];
      A[3] = I[6*i+2];
      A[4] = I[6*i+3];
      A[5] = I[6*i+4];
      A[6] = I[6*i+5];
      triangle(A, J);

      cout<<"B: S = "<<J[0]<<", S' = "<<J[1]<<", S'' = "<<J[3]<<endl;

      d0 += deg0[i];
      d1 += deg1[i];
      d2 += deg2[i];
      

      //cout<<"MU = ("<<I[6*i]<<", "<<I[6*i+1]<<", "<<I[6*i+2]<<")/"<<D[i]<<endl; 
      //cout<<"NU = ("<<I[6*i+3]<<", "<<I[6*i+4]<<", "<<I[6*i+5]<<")/"<<D[i]<<"\n\n";
      //      cout<<J[0]<<J[1]<<J[2]<<J[3]<<J[4]<<J[5]<<endl;

      if(B==J){
	divtemp = square(B,A);
        d0-=divtemp;
        d1-=divtemp;
        d2-=divtemp;
      }
	
      else{
	divtemp = compose(B,J,A);
        d0-=divtemp;
        d1-=divtemp;
        d2-=divtemp;
      }

      cout<<"A*B: S = "<<A[0]<<", S' = "<<A[1]<<", S'' = "<<A[3]<<endl;

      d = A[0];
      ddeg = deg(d);
      mu0 = A[1]*A[2];
      mu1 = A[1];
      mu2 = 0;
      nu0 = A[3]*A[4];
      nu1 = A[3]*A[5];
      nu2 = A[3];
      reduce(half);
      
      cout<<"Basis reduction complete. Reducing the ideal."<<endl;

      /*
	if(zetadeg(mu0, mu1, mu2) >= 0)
	mu0 -= half*trunczeta(mu0, mu1, mu2);
	
	if(zetadeg(nu0, nu1, nu2) >= 0)
	nu0 -= half*trunczeta(nu0, nu1, nu2);
      */
      // The basis may be reduced but the ideal itself may not be.
      // Reduce the ideal.
      reduceideal(half, U);
       
      A[0] = d;
      A[1] = mu0;
      A[2] = mu1;
      A[3] = mu2;
      A[4] = nu0;
      A[5] = nu1;
      A[6] = nu2;

      triangle(A, B);
      cout<<"reduce(A*B): S = "<<B[0]<<", S' = "<<B[1]<<"S'' = "<<B[3]<<endl;
      
      cout<<"At Degree = ("<<d0<<", "<<d1<<", "<<d2<<").\n";

      /*
      if(d0 >=dege1)
	d0%=dege1;
      if(d2 >=dege2)
	d2%=dege2;
      
      cout<<"At Degree = ("<<d0<<", "<<d2<<").\n";
      */
      /*
	cout<<"Reduced?"<<endl;
	cout<<zetadeg(mu0, mu1, mu2)<<" < 0,  "<<zetadeg(nu0, nu1, nu2)<<" <= 0"<<endl;
	cout<<xideg(mu1, mu2)<<" > "<<xideg(nu1, nu2)<<endl;
	cout<<etadeg(mu1, mu2)<<" < 0 <= "<<etadeg(nu1, nu2)<<endl;
	cout<<deg(LeftShift(nu0, shift)+nu1*rho+nu2*omega)-shift-ddeg<<" != 0 "<<endl;
	
	cout<<"Reduced ideal?"<<endl;
	cout<<deg(LeftShift(mu0, shift)+mu1*rho+mu2*omega)- shift-ddeg<<" > 0 "<<endl;
	cout<<"MAX{"<<deg(LeftShift(nu0, shift)+nu1*rho+nu2*omega)- shift-ddeg<<", "<<etadeg(nu1, nu2)<<"} > 0 "<<endl;
      */
      if( search(0, num, 1) == -1 ){
	insert(-l, m, 1);
	m++;
      }
    }
    
    // Square the current ideal.
    else if(step == 6){
      A[0] = d;
      A[1] = mu0;
      A[2] = mu1;
      A[3] = mu2;
      A[4] = nu0;
      A[5] = nu1;
      A[6] = nu2;

      triangle(A, B);
      d0*=2;
      d1*=2;
      d2*=2;

      divtemp = square(B, A);
      d0-=divtemp;
      d1-=divtemp;
      d2-=divtemp;
      
      d = A[0];
      ddeg = deg(d);
      mu0 = A[1]*A[2];
      mu1 = A[1];
      mu2 = 0;
      nu0 = A[3]*A[4];
      nu1 = A[3]*A[5];
      nu2 = A[3];
      reduce(half);

      /*
      if(zetadeg(mu0, mu1, mu2) >= 0)
	mu0 -= half*trunczeta(mu0, mu1, mu2);
      
      if(zetadeg(nu0, nu1, nu2) >= 0)
	nu0 -= half*trunczeta(nu0, nu1, nu2);
      */
      // The basis may be reduced but the ideal itself may not be.
      // Reduce the ideal.
      
      reduceideal(half, U);
      cout<<"At Degree = ("<<d0<<", "<<d1<<", "<<d2<<").\n";
      /*
      if(d0 >=dege1)
	d0%=dege1;
      if(d2 >=dege2)
	d2%=dege2;
      cout<<"At Degree = ("<<d0<<", "<<d2<<").\n";
      */
      /*
	cout<<"Reduced?"<<endl;
	cout<<zetadeg(mu0, mu1, mu2)<<" < 0,  "<<zetadeg(nu0, nu1, nu2)<<" <= 0"<<endl;
	cout<<xideg(mu1, mu2)<<" > "<<xideg(nu1, nu2)<<endl;
	cout<<etadeg(mu1, mu2)<<" < 0 <= "<<etadeg(nu1, nu2)<<endl;
	cout<<deg(LeftShift(nu0, shift)+nu1*rho+nu2*omega)-shift-ddeg<<" != 0 "<<endl;
	
	cout<<"Reduced ideal?"<<endl;
	cout<<deg(LeftShift(mu0, shift)+mu1*rho+mu2*omega)- shift-ddeg<<" > 0 "<<endl;
	cout<<"MAX{"<<deg(LeftShift(nu0, shift)+nu1*rho+nu2*omega)- shift-ddeg<<", "<<etadeg(nu1, nu2)<<"} > 0 "<<endl;
        cout<<deg( LeftShift(nu0, shift) + nu1*rho + nu2*omega)-shift-ddeg<<deg( LeftShift(nu0, shift) + U*nu1*rho + nu2*U*U*omega)-shift-ddeg<<deg( LeftShift(nu0, shift) + U*U*nu1*rho + nu2*U*omega)-shift-ddeg<<endl;
      */

      if( search(0, num, 1) == -1 ){
	insert(-l, m, 1);
	m++;
      }
    }

    else{
      cout<<"ERROR!\n\n";
    }
  }
}
