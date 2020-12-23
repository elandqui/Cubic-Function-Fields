#include <NTL/ZZ_pX.h>
using namespace std;

NTL_CLIENT

// This function will compute the triangle basis of an ideal.
// The input is an ideal A = {1, mu, nu}/d = [d, mu0, mu1, ..., nu2]
// The return value is the vector [s, s', u, s", v, w] where
// I = {1, s'(u + rho), s"(v + w*rho + omega)}/s}

void triangle(vec_ZZ_pX A, vec_ZZ_pX &I){
  ZZ_pX a1, b1, a, b, t, u, w;

  I[0] = A[0];
  XGCD(I[3], a1, b1, A[4], A[6]);
  I[1] = (A[2]*A[6] - A[5]*A[3])/I[3];
  t = ((a1*A[2] + b1*A[5])%I[3])/I[1];
  a = a1 - t*A[6]/I[3];
  b = b1 + t*A[3]/I[3];
  u = I[2] = (A[1]*A[6] - A[4]*A[3])/(I[1]*I[3]);
  I[4] = (a*A[1] + b*A[4])/I[3];
  w = I[5] = (a*A[2] + b*A[5])/I[3];

  I[2] %= (A[0]/I[1]);
  I[5] %= I[1];
  I[4] += u*(I[5]-w);
  I[4] %= (A[0]/I[3]);

  return;
}

// This function finds I = A*B.
// I is in a triangle basis.
// Format: {s, s'(u + rho), s"(v + w rho + omega)}
// [s, s', u, s", v, w]

void compose(vec_ZZ_pX A, vec_ZZ_pX B, vec_ZZ_pX &I){

  ZZ_pX g, ra, rb, sa, sb, u, w;

  if(!IsOne(GCD(A[0], B[0]))){
    cout<<A[0]<<B[0]<<GCD(A[0],B[0])<<endl;
    cout<<"The leading coefficients are not coprime.\n\n";
    return;
  }

  I[0] = A[0]*B[0];  // s_3 = s_1 s_2
  I[1] = A[1]*B[1];  // s_3' = s_1' s_2'
  I[3] = A[3]*B[3];  // s_3" = s_1" s_2"
  
  // Run CRT on the u's, v's, and w's.
  XGCD(g, ra, rb, sa = A[0]/A[1], sb = B[0]/B[1]);
  u = I[2] = (A[2]*rb*sb + B[2]*ra*sa)%(sa*sb);

  XGCD(g, ra, rb, A[1], B[1]);
  w = I[5] = (A[5]*rb*B[1] + B[5]*ra*A[1])%(A[1]*B[1]);

  XGCD(g, ra, rb, sa = A[0]/A[3], sb = B[0]/B[3]);
  I[4] = ((A[4] + A[2]*(I[5]-A[5]))*rb*sb + (B[4] + B[2]*(I[5]-B[5]))*ra*sa)%(sa*sb);

  // "Reduce."

  I[2] %= (I[0]/I[1]);
  I[5] %= I[1];
  I[4] += u*(I[5]-w);
  I[4] %= (I[0]/I[3]);

  return;
}

// This function finds I such that I = A^2.
// The ideal A must have a triangle basis. 
// The ideal I has a triangle basis.

void square(vec_ZZ_pX A, vec_ZZ_pX &I) {
 
  ZZ_pX sg, sh, sg1, sgsh, sonsgh, sonsgh2, umod, g, a, b, sa, sb, temp, temp2;
  ZZ_pX y, z, u, w;
  ZZ_pX H = ZZ_pX(1,1) + ZZ_pX(0,1);
  ZZ_pX G = H + ZZ_pX(4,1); 
  ZZ_pX f;
  ZZ_pX F = G*H*H;

  GCD(sg, A[0], G);
  GCD(sh, A[0], H);
  GCD(sg1, A[1], G);

  sgsh = sg*sh;
  sonsgh = A[0]/sgsh;
  sqr(sonsgh2, sonsgh);
  umod = sonsgh*sg1/A[1];
  
  f = sgsh/(sg1*A[3]);
  //cout<<F<<endl;

  XGCD(g, a, y, umod, 3*sqr(A[2]));
  temp2 = 2*A[5] + H*sqr(A[5]);
  XGCD(g, a, z, sonsgh, temp2);
  if (!IsOne(g)){
    //f = ZZ_pX(0,1);
    A[5] += A[1];
    A[4] += A[2]*A[1];
    XGCD(g, a, z, sonsgh, temp2);
  }

  I[0] = sqr(A[0])/(sgsh);
  I[1] = sqr(A[1])*sg/power(sg1,3);
  I[3] = sh/A[3];

  XGCD(g, a, b, sa = sh*sg1, sb = sqr(umod));
  u = I[2] = (A[2] - y*(power(A[2],3)+F))*a*sa % (sa*sb);
  
  temp = H*power(A[5],3) - G;
  XGCD(g, a, b, sa = sg/sg1, sb = sqr(A[1]/sg1));
  w = I[5] = (A[5] - z*temp)*a*sa % (sa*sb);

  XGCD(g, a, b, sa = sg*A[3], sonsgh2);
  I[4] = (A[4] + I[2]*(I[5]-A[5]) + z*(I[2]*temp + 2*G*H*A[5] - A[4]*(temp2-A[4])))*a*sa % (sa*sonsgh2);

  // "Reduce."
  
  I[2] %= (I[0]/I[1]);
  I[5] %= I[1];
  I[4] += u*(I[5]-w);
  I[4] %= (I[0]/I[3]);

  return;

}

main(){
  vec_ZZ_pX I, A, B;
  ZZ_pX d, mu0, mu1, mu2, nu0, nu1, nu2;

  ZZ_p::init(to_ZZ(2));
  
  I.SetLength(6);
  A.SetLength(7);
  B.SetLength(6);

  A[0] = ZZ_pX(1, 1);
  A[1] = ZZ_pX(0,1) + A[0];
  A[2] = ZZ_pX(2,1) + A[1];
  A[3] = ZZ_pX(0,1);
  A[4] = A[3];
  A[5] = ZZ_pX(3,1) + A[1];
  A[6] = A[1];
  cout<<"A = {"<<A[0]<<", "<<A[1]<<" + "<<A[2]<<"RHO + "<<A[3]<<"OMEGA, "<<A[4]<<" + "<<A[5]<<"RH0 + "<<A[6]<<"OMEGA}\n";
  triangle(A, I);

  cout<<"Triangle basis = {"<<I[0]<<", "<<I[1]<<"("<<I[2]<<" + RHO), "<<I[3]<<"("<<I[4]<<" + "<<I[5]<<"RH0 + OMEGA)}\n\n";
  A.SetLength(6);
  A[0] = ZZ_pX(2, 1);
  A[1] = ZZ_pX(0,1);
  A[2] = A[1] + ZZ_pX(1,1);
  A[3] = A[1];
  A[4] = A[2];


  B[0] = A[2]*(A[2] + ZZ_pX(4,1));
  B[1] = A[2] + ZZ_pX(4,1);
  B[3] = A[2];

  cout<<"A = {"<<A[0]<<", "<<A[1]<<"("<<A[2]<<" + RHO), "<<A[3]<<"("<<A[4]<<" + "<<A[5]<<"RH0 + OMEGA)}\n";
cout<<"B = {"<<B[0]<<", "<<B[1]<<"("<<B[2]<<" + RHO), "<<B[3]<<"("<<B[4]<<" + "<<B[5]<<"RH0 + OMEGA)}\n\n";
  compose(A, B, I);

  cout<<"I = A*B = {"<<I[0]<<", "<<I[1]<<"("<<I[2]<<" + RHO), "<<I[3]<<"("<<I[4]<<" + "<<I[5]<<"RH0 + OMEGA)}\n\n";

  A[0] = ZZ_pX(1,1)*B[0];
  A[1] = A[2] + ZZ_pX(4,1);
  A[3] = A[2];
  A[4] = A[1];

  cout<<"A = {"<<A[0]<<", "<<A[1]<<"("<<A[2]<<" + RHO), "<<A[3]<<"("<<A[4]<<" + "<<A[5]<<"RH0 + OMEGA)}\n";

  square(A, I);
  
cout<<"I = A^2 = {"<<I[0]<<", "<<I[1]<<"("<<I[2]<<" + RHO), "<<I[3]<<"("<<I[4]<<" + "<<I[5]<<"RH0 + OMEGA)}\n\n";
}
