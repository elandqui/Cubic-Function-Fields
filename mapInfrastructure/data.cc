#include "data.h"

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
inline int xideg2(ZZ_pX b, ZZ_pX c, ZZ_p U){
  //cout<<"xi: "<<b*rho + c*omega<<endl;
  return(deg(b*U*U*rho + c*U*omega) - shift-ddeg);
}
// eta_alpha = b*rho - c*omega
inline int etadeg2(ZZ_pX b, ZZ_pX c, ZZ_p U){
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
inline ZZ_pX divxi2(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2, ZZ_p U){
  ZZ_pX ximu, xinu, q;
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
inline ZZ_pX diveta2(ZZ_pX mu0, ZZ_pX mu1, ZZ_pX mu2, ZZ_pX nu0, ZZ_pX nu1, ZZ_pX nu2, ZZ_p U){
  ZZ_pX etamu, etanu, q;
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
inline ZZ_pX trunczeta2(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U){
  ZZ_pX q;
  div(q, LeftShift(2*a, shift) - b*U*U*rho - c*U*omega, LeftShift(d, shift));

  //cout<<"trunczeta2: "<<q<<endl;

  return(q*d);
}
inline ZZ_pX aut(ZZ_pX a, ZZ_pX b, ZZ_pX c, ZZ_p U){
  return(a + RightShift(U*b*rho + U*U*c*omega,shift));
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
	
	mu0 -= lc*nu0; mu1 -= lc*nu1; mu2 -= lc*nu2; 
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
void reduce2(ZZ_p U, ZZ_p half){
  int xs, xt, es, et; 
  ZZ_pX temp0, temp1, temp2, q, etasig, etatau;
  ZZ_p lc, lc1;

  //Step 2 in 4.6
  xs = xideg2(sig1, sig2, U);
  xt = xideg2(tau1, tau2, U);

  if(xs < xt){
    temp0 = sig0; temp1 = sig1; temp2 = sig2;
    sig0 = tau0;  sig1 = tau1; sig2 = tau2;
    tau0 = -temp0; tau1 = -temp1; tau2 = -temp2;
  }
  if(xs == xt){
    es = etadeg2(sig1, sig2, U); et = etadeg2(tau1, tau2, U);
    if(es < et){
      temp0 = sig0; temp1 = sig1; temp2 = sig2;
      sig0 = tau0; sig1 = tau1; sig2 = tau2;
      tau0 = -temp0; tau1 = -temp1; tau2 = -temp2;
    }
  }

  // Step 3.
  es = etadeg2(sig1, sig2, U); et = etadeg2(tau1, tau2, U);
    if(es >= et){
      // Step 3.1
      //cout<<"nu "<<tau1<<tau2<<endl;
       while( deg(LeftShift(sqr(tau1),1)*U*omega - sqr(tau2)*G*U*U*rho) - shift  > (double)discdeg(sig1, sig2, tau1, tau2)/2  ){
	q = divxi2(sig0, sig1, sig2, tau0, tau1, tau2, U);
	temp0 = sig0; temp1 = sig1; temp2 = sig2;
	sig0 = tau0; sig1 = tau1; sig2 = tau2;
	tau0 = q*sig0-temp0; tau1 = q*sig1-temp1; tau2 = q*sig2-temp2;  
      }

      // Step 3.2
      q = divxi2(sig0, sig1, sig2, tau0, tau1, tau2, U);
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

    while(etadeg2(sig1, sig2, U) >= 0){
      q = diveta2(tau0, tau1, tau2, sig0, sig1, sig2, U);

      temp0 = tau0; temp1 = tau1; temp2 = tau2;
      tau0 = sig0; tau1 = sig1; tau2 = sig2;
      sig0*=q; sig1*=q; sig2*=q;
      sig2-=temp2; sig1-=temp1; sig0-=temp0;
    }
    // Step 5.
  
    sig0 -= half*trunczeta2(sig0, sig1, sig2, U);
    tau0 -= half*trunczeta2(tau0, tau1, tau2, U);
}

void quickreg(ZZ_p U, ZZ_p half ){
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
  reduce(half);

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
    if(deg(autnu = aut(nu0, nu1, nu2, U)) == ddeg){
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
    reduce(half);

    //Normalize d;
    d*=inv(LeadCoeff(d));
   
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
  //  cout<<"Step 3 done"<<endl;

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
  //  cout<<"Step 9"<<endl;

  m=0;
  cont = 1;
  // Step 10.
  do{
    // Step 10.1
    reduce2(U, half);

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
      //  cout<<"Step 10"<<endl;

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
    reduce(half);

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
  //cout<<"Step 11"<<endl;

  // Step 12
  for(i=p; i<j; i++){
    // cout<<i<<". "<<N[2*i]<<" "<<N[2*i+1]<<endl;
    e21-=N[2*i];
    e22-=N[2*i+1];
  }
  //cout<<"Step 12"<<endl;

  // Step 13.
  R= abs(e11*e22-e12*e21);

  jj = j-p;

  /*
  cout<<G<<endl;
  cout<<"The preperiod is:    p = "<<p<<endl;
  cout<<"The period is:       l = "<<l<<endl;
  cout<<"The shift is:        j = "<<j-p<<endl;
  cout<<"The period is:       m = "<<m<<endl;
  //cout<<"deg(e1) = "<<e11<<" deg(e1') = "<<e12<<endl;
  //cout<<"deg(e2) = "<<e21<<" deg(e2') = "<<e22<<endl;
  cout<<"The Regulator is:    R = "<<R<<".\n\n";
  */
}

// Stores elements in a hashtable. Returns 0 if there is a duplicate.

int store(int storeR[], int storeP[], int storeL[], int storeM[], int storeJ[], ZZ_pX storeG[], int hashsize){
  int hash = R%hashsize;
  int i=1;

  while(storeR[hash]){
    // If there is a duplicate, don't insert it.
    if((storeR[hash] == R) && (storeP[hash] == p) && (storeL[hash] == l) && (storeM[hash] == m) && (storeJ[hash] == jj))
      return 0;

    // Otherwise, we'll do quadratic probing.
    hash = (hash + i*i)%hashsize;

  }
  storeR[hash] = R;
  storeP[hash] = p;
  storeL[hash] = l;
  storeM[hash] = m;
  storeJ[hash] = jj;
  storeG[hash] = G;

}

void qsort(int storeR[], int storeP[], int storeL[], int storeM[], int storeJ[], ZZ_pX storeG[], int lower, int upper){
  int i, j;
  int temp, pivot;

    if (lower < upper) {
	SWAP(storeR[lower], storeR[(upper+lower)/2]);
	SWAP(storeP[lower], storeP[(upper+lower)/2]);
	SWAP(storeL[lower], storeL[(upper+lower)/2]);
	SWAP(storeM[lower], storeM[(upper+lower)/2]);
	SWAP(storeJ[lower], storeJ[(upper+lower)/2]);
	tempG = storeG[lower]; 
	storeG[lower] = storeG[(upper+lower)/2];
	storeG[(upper+lower)/2] = tempG;
	//SWAP(storeG[lower], storeG[(upper+lower)/2]);

	i = lower;  j = upper + 1;  pivot = storeR[lower];
	while (1) {
	  
	  do i++; while (storeR[i] < pivot);
	  do j--; while (storeR[j] > pivot);
	  if (j < i) break;
	  SWAP(storeR[i], storeR[j]);
	  SWAP(storeP[i], storeP[j]);
	  SWAP(storeL[i], storeL[j]);
	  SWAP(storeM[i], storeM[j]);
	  SWAP(storeJ[i], storeJ[j]);
	  tempG = storeG[i]; 
	  storeG[i] = storeG[j];
	  storeG[j] = tempG;
	  //SWAP(storeG[i], storeG[j]);
	}
	SWAP(storeR[lower], storeR[j]);
	SWAP(storeP[lower], storeP[j]);
	SWAP(storeL[lower], storeL[j]);
	SWAP(storeM[lower], storeM[j]);
	SWAP(storeJ[lower], storeJ[j]);
	tempG = storeG[lower]; 
	storeG[lower] = storeG[j];
	storeG[j] = tempG;
	// SWAP(storeG[lower], storeG[j]);

	qsort (storeR, storeP, storeL, storeM, storeJ, storeG, lower, j - 1);
	qsort (storeR, storeP, storeL, storeM, storeJ, storeG, i, upper);
    }


}

void sort(int storeR[], int storeP[], int storeL[], int storeM[], int storeJ[], ZZ_pX storeG[], int hashsize, int trials){

  int i=0, j=1;

  cout<<"Sorting.\n\n";

  // First we'll group all the elements together.
  while( j < hashsize ){
    if(storeR[i]){
      i++;
      j++;
    }
    else{
      while (!storeR[j]) { j++; }
      storeR[i] = storeR[j];
      storeP[i] = storeP[j];
      storeL[i] = storeL[j];
      storeM[i] = storeM[j];
      storeJ[i] = storeJ[j];
      storeG[i] = storeG[j];
      storeR[j] = 0;
      i++;
      j++;
    }
  }

  cout<<"Done shifting everything. Running Quicksort.\n\n";

  // The largest possible R.
  storeR[trials] = to_int(power(sqrt(to_RR(P)) + 1,6) + 1);
 
  qsort(storeR, storeP, storeL, storeM, storeJ, storeG, 0, trials-1);


  for(i=0; i<trials; i++)
    cout<<storeR[i]<<" "<<storeP[i]<<" "<<storeL[i]<<" "<<storeM[i]<<" "<<storeJ[i]<<" "<<storeG[i]<<endl;
    
  cout<<endl;
}


int main(int argc, char *argv[]){

  int h, i, j, k, n, o;
  int dup; // A duplication flag
  ZZ_pX autnu;
  int coefg[5]; //= {1,1,0,2,1};
  // int badinput = 1;
  int trials;
  ZZ_pX DG; // G'(x)
  int hashsize;

  // Get the prime.
  P = atoi(argv[1]);
  if(P%3 != 1){
    cout<<"The characteristic must be P = 1 (mod 3).\n\n";
    return 0;
  }
  

  // Set the characterisitc.
  ZZ_p::init(to_ZZ(P));

  // Create hashtables
  hashsize = NextPrime(3*P*P*P);
  cout<<"Hashsize = "<<hashsize<<endl;
  int storeR[hashsize];
  int storeP[hashsize];
  int storeL[hashsize];
  int storeM[hashsize];
  int storeJ[hashsize];
  ZZ_pX storeG[hashsize];

  for(i=0; i<hashsize; i++)
    storeR[i] = 0;

  ZZ_p half = inv(to_ZZ_p(2));

  // To test transformations:
  ZZ_p a; 
  ZZ_p mults[P-2];
  ZZ_p halfs[P-2];
  for(i=0; i<(P-2); i++){
    mults[i] = power(to_ZZ_p(i+2), 5);
    halfs[i] = inv(to_ZZ_p(i+2));
  }

  // Size for arrays that store infrastructure ideals.
  maxsize = (int)(8*pow(sqrt((double)P)+1, 4));

  // Number of trials that we run to get data.
  trials = 0;
  

  // Find a primitive cube root of 1
  for(i=1; i<P; i++){
    if( (i*i + i + 1)%P == 0 ) {
      u=i;
      break;
    }
  }
  
  ZZ_p U = to_ZZ_p(u);

  // Find data for all the possible G(x) = x^4 + hx^3 + ix^2 + jx + k.
  SetCoeff(G, 4);
  coefg[4] = 1;
  SetCoeff(DG, 3, 4);
  for(h=0; h<P; h++){
    SetCoeff(G, 3, h);
    coefg[3] = h;
    SetCoeff(DG, 2, 3*h);
    for(i=0;  i<P; i++){
      SetCoeff(G, 2, i);
      coefg[2] = i;
      SetCoeff(DG, 1, 2*i);
      for(j=0; j<P; j++){
        SetCoeff(G, 1, j);
	coefg[1] = j;
	SetCoeff(DG, 0, j);
        for(k=1; k<P; k++){
          SetCoeff(G, 0, k);
	  coefg[0] = k;

          dup = 0;
          // Since G(ax)H^2(ax) = G(x)H^2(x), we get duplicate info
          // for all P-1 G(x), G(2x), G(3x), ..., G(-x).
          // Only choose those G whose coefficients are the smallest
          // under this transformation.
          for(n=2; (n<(P-1)) && !dup ; n++){
            a = mults[n-2];
	    for(o=3; o>=0 ; o--){
	      if(coefg[o] == 0){
	        a*=halfs[n-2];
	        continue;
	      }
	      if(rep(coeff(G, o)) > rep(a*coeff(G, o))){
		dup = 1;
		break;
	      }
	      if(rep(coeff(G, o)) < rep(a*coeff(G,o))){
		break;
	      }
	      a*=halfs[n-2];
	    }
          }

	  // We only want G if there are no repeated roots.
	  // This is true iff gcd(G, G') = 1.
	  if(IsOne(GCD(DG,G)) && !dup){

            cout<<trials<<" "<<G<<endl;

	    // Run the tests on this polynomial.
	    // Find a good approximation for the basis elements.
	    basis(shift);
	    
	    // The coefficients of the basis elements.
	    mu0 = ZZ_pX(0,0), mu1 = ZZ_pX(0,1), mu2 = ZZ_pX(0,0); 
	    nu0 = ZZ_pX(0,0), nu1 = ZZ_pX(0,0), nu2 = ZZ_pX(0,1);
	    
	    // Set the denominator
	    d = to_ZZ_pX(1);
	    ddeg = 0;

	    p=0; l=0; m=0;

	    // Run the reduction.
	    reduce(half);
	    
	    // Find the regulator.
	    quickreg(U, half);
	    
            if(trials + 5 > hashsize)
              cout<<"Getting close! "<<trials<<G<<endl;
            
	    // Store the data.
	    if(store(storeR, storeP, storeL, storeM, storeJ, storeG, hashsize)){
	      trials++;
              Rsum+=R;
              psum+=p;
              lsum+=l;
              jsum+=jj;
              msum+=m;
              if(!p)
                zerop++;
	    }

	  }
	  // else continue.
        }
      }
    }
  }

  cout<<"Done now. Found "<<trials<<" polynomials.\n\n";

  // Sort the hash table by R and print it out.
  sort(storeR, storeP, storeL, storeM, storeJ, storeG, hashsize, trials);

  cout<<"Total runs: "<<trials<<endl;
  cout<<"Times p = 0: "<<zerop<<endl;
  cout<<"Average p = "<<psum/trials<<endl;
  cout<<"Average l = "<<lsum/trials<<endl;
  cout<<"Average m = "<<msum/trials<<endl;
  cout<<"Average j = "<<jsum/trials<<endl;
  cout<<"Average R = "<<Rsum/trials<<endl;
  //  cout<<"Total p = "<<psum<<endl;
  //cout<<"Total l = "<<lsum<<endl;
  //cout<<"Total m = "<<msum<<endl;
  //cout<<"Total j = "<<jsum<<endl;
  //cout<<"Total R = "<<Rsum<<endl;
}
