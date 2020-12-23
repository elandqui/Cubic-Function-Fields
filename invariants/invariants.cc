/*
  invariants.cc 

  Contains global variables for use by any program computing with 
  cubic function fields.

  To make the library libinvariants, put invariants.cc, invariants.h, and
  makefile into a directory, Invariants. Make sure the path to NTL and GMP 
  are set correctly in the makefile and run make.

  Written by Eric Landquist.

 */

#include "invariants.h"

NTL_CLIENT

ZZ q;        // The characteristic.

ZZ_pX G;     // The defining polynomials.
ZZ_pX H;     // The defining polynomials.
ZZ_pX f;     // f = GH^2. The function field is F_q(C), C: Y^3 = f.
int genus=3; // The genus of the function field.
int rank=0;  // The unit rank of F_q[C].

int NUMPRIMES=320;

long PRIMES[320] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823,  827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129};

// The following constants are based on experimental evidence.
// The constants are the average values of |E-h|/U, where:
// h is the divisor class number of F_q(C).
// E is the estimate of h
// U is the upper bound on the error
// E and U are found via the method of Scheidler and Stein, 2008:
// "Approximating Euler products and class number computation 
//                in algebraic function fields"
//
// The else is a pure guess.

double alphahat(){ 
  if (genus == 3)
    return 0.27187490;
  else if (genus == 4)
    return 0.19186318;
  else if (genus == 5)
    return 0.19190607; 
  else if (genus == 6)
    return 0.15975657; 
  else if (genus == 7)
    return 0.12602172;
  else
    return 1.0/(2.0*((double)genus-1.0)); 

}

// Ratio between a giant step and inverse and a baby step in rank 0.

double tau1(){ 
  if (genus == 3)
    return 1.38; 
  else if (genus == 4)
    return 1.51; 
  else if (genus == 5)
    return 1.393; 
  else if (genus == 6)
    return 1.615;
  else if (genus == 7)
    return 1.54;
  else 
    return 1.615;
}

// Ratio between a giant step and inverse and a baby step in rank 1.

double tau2(){ 
  if (genus == 3){
    if(deg(G) == 1)
      return 4.47472;
    else // if deg(G) == 4.
      return 4.20202;
  }
  else if(genus == 4){
    if(deg(G) == 3)
      return 6.11146;
    else // if deg(G) == 6.
      return 5.49073;
  }
  else if(genus == 5){
    if(deg(G) == 2)
      return 8.47669;
    else // if deg(G) == 5.
      return 8.09893;
  }
  else if(genus == 6){
    if(deg(G) == 1)
      return 9.69945;
    else if(deg(G) == 4)
      return 9.33196;
    else // if deg(G) == 7.
      return 9.02126;
  }
  else if(genus == 7){
    if(deg(G) == 3)
      return 11.90331;
    else if(deg(G) == 6)
      return 11.56116;
    else // if deg(G) == 9.
      return 11.26627;
  }
  else
    return ((double)genus)*1.5;
} 

// Ratio between a giant step and a baby step in rank 1.

double tau3(){
  if(genus == 2){
    return 2.96977;
  }
  else if(genus == 3){
    if(deg(G) == 1)
      return 3.17201;
    else // if deg(G) == 4.
      return 2.92374; 
  }
  else if(genus == 4){
    if(deg(G) == 3)
      return 4.11812;
    else // if deg(G) == 6.
      return 3.87316;
  }
  else if(genus == 5){
    if(deg(G) == 2)
      return 5.57365;
    else // if deg(G) == 5.
      return 5.29813;
  }
  else if(genus == 6){
    if(deg(G) == 1)
      return 6.64731;
    else if(deg(G) == 4)
      return 6.10144;
    else // if deg(G) == 7.
      return 5.86166;
  }
  else if(genus == 7){
    if(deg(G) == 3)
      return 7.85012;
    else if(deg(G) == 6)
      return 7.72477;
    else // if deg(G) == 9.
      return 7.50799;
  }
  else
    return (double)genus;
}

// Ratio between a giant step and inverse and a baby step in rank 2.
// In unit rank 2, using an inverse is slower than using a second
// giant step because the inverse requires additional baby steps, 
// including basis reductions, to get back on track.
// These extra steps prove to be too costly.
double tau4(){ 
  return (2.0*tau5());
} 

// Ratio between a giant step and a baby step in rank 2 (estimated).

double tau5(){
  if(genus == 2){
    return 3.04839;  
  }
  else if(genus == 3){
    if(deg(G) == 1)
      return 3.28523;  
    else // if deg(G) == 4.
      return 3.02410;  
  }
  else if(genus == 4){
    if(deg(G) == 3)
      return 4.42018; 
    else // if deg(G) == 6.
      return 4.03846; 
  }
  else if(genus == 5){
    if(deg(G) == 2)
      return 5.93367; 
    else // if deg(G) == 5.
      return 5.61416; 
  }
  else if(genus == 6){
    if(deg(G) == 1)
      return 7.10509; 
    else if(deg(G) == 4)
      return 6.38660;
    else // if deg(G) == 7.
      return 5.96440;
  }
  else if(genus == 7){
    if(deg(G) == 3)
      return 8.58138; 
    else if(deg(G) == 6)
      return 8.21655; 
    else // if deg(G) == 9.
      return 7.87264; 
  }
  else
    return (double)genus;
}

// For use with the kangaroo method in unit rank 1:
// Determines the "headwind" experienced by kangaroos via
// reduction after a giant step with probability 1-O(1/q).
// That is, if A = B*C, head0 = dist(B) + dist(C) - dist(A)
// in almost all cases.
// Input: genus - the genus of F_q(C).
// Output: head0 - the headwind as described above.
void headwind1(long &head0){
  if(genus%3 == 1)
    head0 = (genus+2)/3;
  else
    head0 = genus/3;
}

// For use with the kangaroo method in unit rank 2:
// Determines the "headwind" experienced by kangaroos via
// reduction after a giant step with probability 1-O(1/q).
// That is, if A = B*C, then 
// head0 = dist_0(B) + dist_0(C) - dist_0(A)
// head1 = dist_1(B) + dist_1(C) - dist_1(A)
// head2 = dist_2(B) + dist_2(C) - dist_2(A)
// in almost all cases.
// Input: genus - the genus of F_q(C).
// Output: head0,1,2 - the headwind as described above.
void headwind2(long &head0, long &head1, long &head2){
  if(genus%3 == 1)
    head0 = (genus+2)/3;
  else
    head0 = genus/3;
  if(genus%3 == 2)
    head1 = head2 = (genus+1)/3;
  else
    head1 = head2 = genus/3;
}

// For use with the BSGS method in unit rank 2:
// Determines the 2-distance "headwind" experienced by giant 
// steps via reduction after an inverse with probability 1-O(1/q).
// That is, if A = Inverse(B), then 
// headwind = dist_2(B) - dist_0(A)
// in almost all cases.
// Input: genus - the genus of F_q(C).
// Output: the headwind as described above.
long inverseheadwind2(){
  if(genus%3 == 0)
    return (2*genus/3);
  else
    return ((2*genus + 1)/3);
}

// Used to optimize below() in infrastructure.cc.
// For optimizing the computation of the infrastructure divisor
// whose distance is below a given integer.

int belowsteps(){
  if (genus == 3){
    if(deg(G) == 1)
      return 5;
    else
      return 4;
  }
  else if(genus == 4){
    return 6;
  }
  else if(genus == 5){
    if(deg(G) == 2)
      return 9;
    else // if deg(G) == 5.
      return 8;
  }
  else if(genus == 6){
    if(deg(G) == 1)
      return 10;
    else if(deg(G) == 4)
      return 10;
    else // if deg(G) == 7.
      return 9;
  }
  else if(genus == 7){
    if(deg(G) == 3)
      return 12;
    else if(deg(G) == 6)
      return 11;
    else // if deg(G) == 9.
      return 11;
  }
  else{
    if(deg(H) > deg(G))
      return (int)((37.2850261249*((double)genus) + 117.0)/37.4299477502 + 0.5);
    else
      return (int)((37.2850261249*((double)genus) + 69.0)/37.4299477502 + 0.5);
  }
}
