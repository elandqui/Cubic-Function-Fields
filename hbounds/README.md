# libhbounds library

  This library contains utilities for computing an estimate, E, of the divisor class number, h, of the cubic function field K = Fq(C), C: Y^3 = f = GH^2, where char(K) > 3.

   **hbounds.cc**

   This file contains utilities for computing an estimate, E, of 
   the divisor class number, h, of the cubic function field K = Fq(C), 
   C: Y^3 = f = GH^2, and an upper bound, U, on the error |h-E|. 
   For the functions, we will use results and algorithms from:

   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)
   [SS7] R. Scheidler and A. Stein, "Class number approximation in 
                      cubic function fields"
   [SS8] R. Scheidler and A. Stein, "Approximating Euler products 
                      and class number computation in algebraic 
                      function fields"
   
   How to create the library:
   --------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download hbounds.cc, hbounds.h and makefile 
      into ../Hbounds.
   Compile libhbounds via make in ../Hbounds.

   Make sure that the paths to the NTL libraries are
   specified correctly in the makefile.
   
   Known Problems: 
   ---------------
   none
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   3.a. Possibly apply sieving in Sv1 and Sv2 for lambda > 2.
   3.b. rec3() rec32(), and isCube() can probably be improved via reciprocity.

   Author: 
   Eric Landquist
