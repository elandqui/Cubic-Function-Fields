#include "infmap.h"

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

inline ZZ_pX divxi(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2){
  ZZ_pX ximu, xinu, q;
  ximu = mu1*rho + mu2*omega;
  xinu = nu1*rho + nu2*omega;

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
  ZZ_pX q;
  
  div(q, LeftShift(2*a, shift) - b*rho - c*omega, LeftShift(d, shift));

  return(q*d);
}
inline ZZ_pX trunczeta2(ZZ_pX a, ZZ_pX b, ZZ_pX c){
  ZZ_pX q;
  ZZ_p U = to_ZZ_p(u);
  div(q, LeftShift(2*a, shift) - b*U*U*rho - c*U*omega, LeftShift(d, shift));

  //cout<<"trunczeta2: "<<q<<endl;

  return(q*d);
}
inline ZZ_pX aut(ZZ_pX a, ZZ_pX b, ZZ_pX c){
  ZZ_p U = to_ZZ_p(u);
  return(a + RightShift(U*b*rho + U*U*c*omega,shift));
}

void reduce(){
  int xm, xn, em, en; 
  ZZ_pX temp0, temp1, temp2, q, etamu, etanu;
  ZZ_p lc, lc1;
  ZZ_p half = inv(to_ZZ_p(2));
   
  //Step 2 in 4.6
  xm = xideg(mu1, mu2);
  xn = xideg(nu1, nu2);

  if(xm < xn){
    temp0 = mu0; temp1 = mu1; temp2 = mu2;
    mu0 = nu0;  mu1 = nu1; mu2 = nu2;
    nu0 = -temp0; nu1 = -temp1; nu2 = -temp2;
  }
  if(xm == xn){
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
	
	mu0 -= lc1*nu0; mu1 -= lc1*nu1; mu2 -= lc1*nu2; 
      }
    }  

    //Step 4.
    /*while(etadeg(nu1, nu2) < 0){
      q = divxi(nu0, nu1, nu2, mu0, mu1, mu2);

      temp0 = mu0; temp1 = mu1; temp2 = mu2;
      mu0 = nu0; mu1 = nu1; mu2 = nu2;
      nu0*=q; nu1*=q; nu2*=q;
      nu2-=temp2; nu1-=temp1; nu0-=temp0;
      cout<<"Step 4a: "<<endl;
	cout<<"mu: {" << mu0<<", "<<mu1<<", "<<mu2<<"}"<<endl;
	cout<<"nu: {" << nu0<<", "<<nu1<<", "<<nu2<<"}"<<endl;
	}*/

    while(etadeg(mu1, mu2) >= 0){
      q = diveta(nu0, nu1, nu2, mu0, mu1, mu2);

      temp0 = nu0; temp1 = nu1; temp2 = nu2;
      nu0 = mu0; nu1 = mu1; nu2 = mu2;
      mu0*=q; mu1*=q; mu2*=q;
      mu2-=temp2; mu1-=temp1; mu0-=temp0;
    }  
 
    // Step 5.
  
    mu0 -= half*trunczeta(mu0, mu1, mu2);
    nu0 -= half*trunczeta(nu0, nu1, nu2);

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

// Computes a 2-reduced basis of {1, mu, nu}.
void reduce2(){
  int xs, xt, es, et; 
  ZZ_pX temp0, temp1, temp2, q, etasig, etatau;
  ZZ_p lc, lc1;
  ZZ_p half = inv(to_ZZ_p(2));
  ZZ_p U = to_ZZ_p(u);

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

    while(etadeg2(sig1, sig2) >= 0){
      q = diveta2(tau0, tau1, tau2, sig0, sig1, sig2);

      temp0 = tau0; temp1 = tau1; temp2 = tau2;
      tau0 = sig0; tau1 = sig1; tau2 = sig2;
      sig0*=q; sig1*=q; sig2*=q;
      sig2-=temp2; sig1-=temp1; sig0-=temp0;
    }
    // Step 5.
  
    sig0 -= half*trunczeta2(sig0, sig1, sig2);
    tau0 -= half*trunczeta2(tau0, tau1, tau2);
}

void babystep1(){
  ZZ_pX a0, a1, a2; // neighbor
  ZZ_pX o0, o1, o2, i0, i1, i2; // temp variables
  int i;
  ZZ_pX dold;  //old denominator
  ZZ_pX autnu;
  ZZ_pX g;            // gcd of numerator and denominator

  // Step 2 is reducing {1, rho, omega}.
  reduce();
  // Step 3.
  // Step 3.1
  
  // Step 3.2
  if(deg(autnu = aut(nu0, nu1, nu2)) == ddeg){
    a0 = nu0-d*ZZ_p(LeadCoeff(autnu));
    a1 = nu1; a2=nu2;
    
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
  reduce();
  
  //Normalize d;
  d*=inv(LeadCoeff(d));
  cout<<"MU_{("<<l<<", "<<m<<")} = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
  cout<<"NU_{("<<l<<", "<<m<<")} = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
  
  
}


void babystep2(){
  int cont = 1;
  int j;
  ZZ_pX dold;
  ZZ_pX o0, o1, o2, i0, i1, i2; // temp variables
  ZZ_pX b0, b1, b2; // beta
  ZZ_pX tau;
  ZZ_pX g; // gcd
  //p=0;

  // Step 2.1
  reduce2();
  
  // Step 2.2
  if(deg(tau = LeftShift(tau0,shift)+tau1*rho+tau2*omega) == ddeg+shift){
    b0 = tau0 - LeadCoeff( tau )*d; b1 = tau1; b2 = tau2;
    
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
  
  reduce();
  
}



int main(){

  int i, j;
  ZZ_pX autnu;
  int coefg[5]; //= {1,1,0,2,1};
  int badinput = 1;
  long ilen, dlen; // clen;
  int num=0;
  int cont = 1;
  int step;
  int match;

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

  // Find a primitive cube root of 1
  for(i=1; i<P; i++){
    if( (i*i + i + 1)%P == 0 ) {
      u=i;
      break;
    }
  }

  // The coefficients of the basis elements.
  t0 = to_ZZ_pX(1);
  mu0 = ZZ_pX(0,0), mu1 = ZZ_pX(0,1), mu2 = ZZ_pX(0,0); 
  nu0 = ZZ_pX(0,0), nu1 = ZZ_pX(0,0), nu2 = ZZ_pX(0,1);

  // Set the denominator
  d = to_ZZ_pX(1);

  // Initialize rho, omega, and G.
  // The coefficients of rho, omega, and G.
  //int coefr[8] = {6,3,3,6,1,5,3,1};
  //int coefo[9] = {4,2,1,3,1,4,5,6,1};
  
  /*for(i = 0; i < 3+shift; i++)
    rho += ZZ_pX(i, coefr[i]);
  for(i = 0; i < 4+shift; i++)
  omega += ZZ_pX(i, coefo[i]); */


  for(i=0; i<4;i++)
    G += ZZ_pX(i,coefg[i]);
  G+=ZZ_pX(4, 1);

  //cout<<"RHO: "<<rho<<endl;
  //cout<<"OMEGA: "<<omega<<endl;
  
  // Find a good approximation for the basis elements.
  basis(5);

  // Run the reduction.
  reduce();

  cout<<"MU_{("<<l<<", "<<m<<")} = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
  cout<<"NU_{("<<l<<", "<<m<<")} = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
 
  l=0;
  cont=1;
  //Get the initial period.
  do{
    //C[3*l] = t0; C[3*l+1]=t1; C[3*l+2]=t2;
    I[6*l] = mu0; I[6*l+1]=mu1; I[6*l+2]=mu2;
    I[6*l+3] = nu0; I[6*l+4]=nu1; I[6*l+5]=nu2;
    D[l] = d;
    ddeg=deg(d);
    coordl[l] = l;
    coordm[l] = 0;
    l++;
    num++;
    babystep1();

    // See if {1, mu, nu} is in I.
    // It's easier to first check for the denominator.
    for(p=0; p< l; p++){
      if(D[p]==d){
	// Now check for mu and nu.
	if((I[6*p] == mu0) && (I[6*p+1] == mu1) && (I[6*p+2] == mu2) && (I[6*p+3] == nu0) && (I[6*p+4] == nu1) && (I[6*p+5] == nu2)){
	  cont=0;
	  break;
	}
      }
    }

  } while(cont);
  
  // Check for the preperiod.
  if(p)
    l-=p;

  cout<<"PREPERIOD p = "<<p<<endl;
  cout<<"PERIOD #1 l = "<<l<<"\n\n"; 

  // Get the secondary periods
  for(i=p; i<p+l; i++){
    m=0;
    //phi0 = C[3*i]; phi1 = C[3*i+1]; phi2 = C[3*i+2];
    sig0 = I[6*i]; sig1=I[6*i+1]; sig2=I[6*i+2];
    tau0=I[6*i+3]; tau1=I[6*i+4]; tau2 = I[6*i+5];
    d = D[i];
    ddeg=deg(d);

    cont=1;

    while(cont){
      m++;
      babystep2();

      for(j=0; j< num; j++){
        if(D[j]==d){
          // Now check for mu and nu.
          if((I[6*j] == mu0) && (I[6*j+1] == mu1) && (I[6*j+2] == mu2) && (I[6*j+3] == nu0) && (I[6*j+4] == nu1) && (I[6*j+5] == nu2)){
            cout<<"PERIOD #2 m_"<<i<<" = "<<m<<" AT "<<endl; 
            cout<<"MATCH = ("<<coordl[j]<<", "<<coordm[j]<<") = \n(";
            cout<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
            cout<<" ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
            cont=0;
            break;
          }
        }
      }
      if(cont){
        cout<<"MU_{("<<i<<", "<<m<<")} = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
        cout<<"NU_{("<<i<<", "<<m<<")} = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
        //C[3*num] = phi0; C[3*num+1]=phi1; C[3*num+2]=phi2;
        I[6*num] = mu0; I[6*num+1]=mu1; I[6*num+2]=mu2;
        I[6*num+3] = nu0; I[6*num+4]=nu1; I[6*num+5]=nu2;
        D[num] = d;
        ddeg=deg(d);
        coordl[num] = i;
        coordm[num] = m;
        num++;
      }
    }
  }
  
  cout<<"TOTAL SO FAR = "<<num<<endl;

  //t0 = C[0]; t1 = C[1]; t2 = C[2];
  mu0 = I[0]; mu1=I[1]; mu2=I[2];
  nu0=I[3]; nu1=I[4]; nu2 = I[5];
  d = D[0];
  ddeg=deg(d);

  cont = 1;
  i=0; m=0;
  while (cont){
     match=0;
    cout<<"3 to quit, 1 to restart, Step 0 or 2? ";
    cin>>step;
    
    if(step == 3){
      cont = 0;
      break;
    }
    else if( step == 1){
      //t0 = C[0]; t1 = C[1]; t2 = C[2];
      mu0 = I[0]; mu1=I[1]; mu2=I[2];
      nu0=I[3]; nu1=I[4]; nu2 = I[5];
      d = D[0];
      ddeg=deg(d);
      l=0;
      m=0;
    }
    else if (step == 0){
      l++;
      babystep1();
      for(i=0; i< num; i++){
        if(D[i]==d){
          // Now check for mu and nu.
          if((I[6*i] == mu0) && (I[6*i+1] == mu1) && (I[6*i+2] == mu2) && (I[6*i+3] == nu0) && (I[6*i+4] == nu1) && (I[6*i+5] == nu2)){
            cout<<"MATCH = ("<<coordl[i]<<", "<<coordm[i]<<") = \n("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
            cout<<" ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
            match = 1;
            break;
          }
        }
      }
      if(!match){
        cout<<"MU_{("<<l<<", "<<m<<")} = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
        cout<<"NU_{("<<l<<", "<<m<<")} = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
        //C[3*num] = phi0; C[3*num+1]=phi1; C[3*num+2]=phi2;
        I[6*num] = mu0; I[6*num+1]=mu1; I[6*num+2]=mu2;
        I[6*num+3] = nu0; I[6*num+4]=nu1; I[6*num+5]=nu2;
        D[num] = d;
        num++;
      }

    }
    else if (step == 2){
      phi0=t0; phi1=t1; phi2=t2;
      sig0=mu0; sig1=mu1; sig2=mu2;
      tau0=nu0; tau1=nu1; tau2=nu2;
      m++;
      babystep2();
for(i=0; i< num; i++){
        if(D[i]==d){
          // Now check for mu and nu.
          if((I[6*i] == mu0) && (I[6*i+1] == mu1) && (I[6*i+2] == mu2) && (I[6*i+3] == nu0) && (I[6*i+4] == nu1) && (I[6*i+5] == nu2)){
            cout<<"MATCH = ("<<coordl[i]<<", "<<coordm[i]<<") = \n("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
            cout<<" ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
            match = 1;
            break;
          }
        }
      }
      if(!match){
        cout<<"MU_{("<<l<<", "<<m<<")} = ("<<mu0<<", "<<mu1<<", "<<mu2<<")/"<<d<<endl; 
        cout<<"NU_{("<<l<<", "<<m<<")} = ("<<nu0<<", "<<nu1<<", "<<nu2<<")/"<<d<<"\n\n";
        //C[3*num] = phi0; C[3*num+1]=phi1; C[3*num+2]=phi2;
        I[6*num] = mu0; I[6*num+1]=mu1; I[6*num+2]=mu2;
        I[6*num+3] = nu0; I[6*num+4]=nu1; I[6*num+5]=nu2;
        D[num] = d;
        num++;
      }

    }
    else{
      cout<<"ERROR!\n\n";
    }
  }
}
