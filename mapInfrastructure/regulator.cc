#include "regulator.h"

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

// Computes a 2-reduced basis of {1, sigma, tau}.
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

void quickreg(){
  int N[maxsize];  // List of degrees. 
  int cont = 1;
  int ilen=0, dlen = 0, nlen=0;
  int e11=0, e12=0, e21=0, e22=0;
  int i, j;
  ZZ_pX a0, a1, a2; // neighbor
  ZZ_pX o0, o1, o2, i0, i1, i2; // temp variables
  ZZ_pX dold;
  ZZ_pX b0, b1, b2; // beta
  ZZ_pX tau, autnu;
  ZZ_pX g; // gcd

  // The coefficients of the basis elements.
  mu0 = ZZ_pX(0,0), mu1 = ZZ_pX(0,1), mu2 = ZZ_pX(0,0); 
  nu0 = ZZ_pX(0,0), nu1 = ZZ_pX(0,0), nu2 = ZZ_pX(0,1);

  // Set the denominator
  d = to_ZZ_pX(1);
  ddeg = deg(d);

  // Step 2
  reduce();

  //Step 3.
  do{
    // Step 3.1
   
    //cout<<dlen<<". mu: {" << mu0<<", "<<mu1<<", "<<mu2<<"}/"<<d<<endl;
    //cout<<"nu: {" << nu0<<", "<<nu1<<", "<<nu2<<"}/"<<d<<endl;

    ilen+=6;
    dlen++;
    I.SetLength(ilen);
    D.SetLength(dlen);

    I[ilen-6] = mu0; I[ilen-5]=mu1; I[ilen-4]=mu2;
    I[ilen-3] = nu0; I[ilen-2]=nu1; I[ilen-1]=nu2;
    D[dlen-1] = d;

    // Step 3.2
    if(deg(autnu = aut(nu0, nu1, nu2)) == ddeg){
      a0 = nu0-d*ZZ_p(LeadCoeff(autnu));
      a1 = nu1; a2=nu2;

      // Step 3.3
      N[nlen] = deg(LeftShift(a0, shift)+a1*rho+a2*omega)-ddeg-shift;
      nlen++;
      N[nlen] = deg(LeftShift(a0, shift)+a1*u*rho+a2*u*u*omega)-ddeg-shift;
      nlen++;

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
      N[nlen] = deg(LeftShift(mu0, shift)+mu1*rho+mu2*omega)-ddeg-shift;
      nlen++;
      N[nlen] = deg(LeftShift(mu0, shift)+mu1*u*rho+mu2*u*u*omega)-ddeg-shift;
      nlen++;

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

    // Step 3.5 Reduce the ideal {1, mu, nu}.
    reduce();

    //Normalize d;
    d*=inv(LeadCoeff(d));
   
    if(!(dlen%10000))
      cout<<"Completed "<<dlen<<" steps.\n";

    // See if {1, mu, nu} is in I.
    // It's easier to first check for the denominator.
    for(p=0; p< dlen-1; p++){
      if(D[p]==d){
	// Now check for mu and nu.
	if((I[6*p] == mu0) && (I[6*p+1] == mu1) && (I[6*p+2] == mu2) && (I[6*p+3] == nu0) && (I[6*p+4] == nu1) && (I[6*p+5] == nu2)){
          
	  cont=0;
	  break;
	}
      }
    }
  }while(cont);

  cout<<"Match found: p = "<<p<<" and l = "<<l<<endl;

  // Step 4: p is set.
  // Step 5: just ignore the first p elements
  // Step 6
  l = dlen-p;
  // Step 7: done
  // Step 8 

  for(i=p; i<p+l; i++){
    //cout<<i<<". "<<N[2*i]<<" "<<N[2*i+1]<<endl;
    e11+=N[2*i];
    e12+=N[2*i+1];
  }

  // Step 9
  sig0 = I[6*p+0], sig1 = I[6*p+1], sig2 = I[6*p+2];
  tau0=I[6*p+3], tau1 = I[6*p+4], tau2 = I[6*p+5];
  d = D[p];
  ddeg = deg(d);

  m=0;
  cont = 1;
  // Step 10.
  do{
    // Step 10.1
    reduce2();

    // Step 10.2
    if(deg(tau = LeftShift(tau0,shift)+tau1*rho+tau2*omega) == ddeg+shift){
      b0 = tau0 - LeadCoeff( tau )*d; b1 = tau1; b2 = tau2;

      //Step 10.3
      e21 += deg(LeftShift(b0, shift)+b1*rho + b2*omega)-ddeg-shift;
      e22 += deg(LeftShift(b0, shift)+b1*u*rho + b2*u*u*omega)-ddeg-shift;

      //Step 10.4
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
      //Step 10.3
      e21 += deg(LeftShift(sig0, shift)+sig1*rho + sig2*omega)-ddeg-shift;
      e22 += deg(LeftShift(sig0, shift)+sig1*u*rho + sig2*u*u*omega)-ddeg-shift;

      //Step 10.4
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

    //Step 10.5
    mu0=sig0; mu1=sig1; mu2=sig2;
    nu0=tau0; nu1=tau1; nu2=tau2;
    reduce();

    //Step 10.6
    m++;

    //cout<<"sig: {" <<sig0<<", "<<sig1<<", "<<sig2<<"}/"<<d<<endl;
    //cout<<"tau: {" <<tau0<<", "<<tau1<<", "<<tau2<<"}/"<<d<<endl;
 
    // See if {1, mu, nu} is in I.
    // It's easier to first check for the denominator.
    for(j=p; j < p+l; j++){
      if(D[j]==d){
	// Now check for mu and nu.
	if((I[6*j] == mu0) && (I[6*j+1] == mu1) && (I[6*j+2] == mu2) && (I[6*j+3] == nu0) && (I[6*j+4] == nu1) && (I[6*j+5] == nu2)){
	  cont=0;
	  break;
	}
      }
    }

  }while(cont);

  // Step 11:

  // Step 12
  for(i=p; i<j; i++){
    // cout<<i<<". "<<N[2*i]<<" "<<N[2*i+1]<<endl;
    e21-=N[2*i];
    e22-=N[2*i+1];
  }

  // Step 13.
  R= abs(e11*e22-e12*e21);

  cout<<"\n";
  cout<<"The preperiod is:    p = "<<p<<endl;
  cout<<"The period is:       l = "<<l<<endl;
  cout<<"The shift is:        j = "<<j<<endl;
  cout<<"The period is:       m = "<<m<<endl;
  cout<<"deg(e1) = "<<e11<<" deg(e1') = "<<e12<<endl;
  cout<<"deg(e2) = "<<e21<<" deg(e2') = "<<e22<<endl;
  

  cout<<"The Regulator is:    R = "<<R<<".\n\n";
}

int main(){

  int i;
  ZZ_pX autnu;
  int coefg[5]; //= {1,1,0,2,1};
  int badinput = 1;
  time_t set, uniti, unitii, reg;
  char yn;
  ZZ seed = to_ZZ( time(NULL) );
  SetSeed( seed );


  while(badinput){
    cout<<"\n";
    cout<<"Characteristic of the base field, P: ";
    cin>>P;

    if(P%3 == 1)
      badinput = 0;
    else{
      cout<<"\n";
      cout<<"The characteristic must be P = 1 (mod 3).\n";
    }
  }

  // Set the characterisitc.
  ZZ_p::init(to_ZZ(P));

  maxsize = (int)(8*pow(sqrt((double)P)+1, 4));

  cout<<"\n";
  cout<<"Generate a random prime polynomial G(x)? y/n: ";
  cin>>yn;

  if((yn == 'y') || (yn == 'Y')){
    G = BuildIrred_ZZ_pX(4);
    for(i=0; i<4; i++){
      coefg[i] = to_int(rep(coeff(G, i)));
    }
    cout<<"\n";
    cout<<"Working with the polynomial \n"<<"x^4 + ";
    for(i=3; i>1; i--){
      if(coefg[i]) {
	if(coefg[i] > 1)
	  cout<<coefg[i];
      	cout<<"x^"<<i<<" + ";
      }
    }
    if(coefg[1]){
      if(coefg[1] > 1)
	cout<<coefg[i];
      	cout<<"x + ";
    }
    cout<<coefg[0]<<".\n";
  }
  else{
    cout<<"\n";
    cout<<"Coefficients of G = g_0 + g_1*x + g_2*x^2 + g_3*x^3 + x^4:\n";
    for( i=0; i<4; i++){
      cout<<"g_"<<i<<": ";
      cin>>coefg[i];
    }
    for(i=0; i<4;i++)
      G += ZZ_pX(i,coefg[i]);
    G+=ZZ_pX(4, 1);
    
  }

  // Find a primitive cube root of 1
  for(i=1; i<P; i++){
    if( (i*i + i + 1)%P == 0 ) {
      u=i;
      break;
    }
  }

  // The coefficients of the basis elements.
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

  //cout<<"RHO: "<<rho<<endl;
  //cout<<"OMEGA: "<<omega<<endl;
  
  // Find a good approximation for the basis elements.
  basis(shift);

  set = time(NULL);

  // Run the reduction.
  reduce();
  
  quickreg();

  reg = time(NULL);

  cout<<"R computed in "<<reg-set<<" seconds.\n\n";

}
