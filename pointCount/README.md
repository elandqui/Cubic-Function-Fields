# Point Count

**bsgstiming.cc**

   This file contains functions to compute timing data for the
   ratio between a baby step and giant step (possibly plus an inverse).
   Timings will involve 10^6 steps each.
   For the functions, we will use results and algorithms from:

   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)
   
   
   How to create the executable:
   -----------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download cubic_ideal.cpp, cubic_ideal.h and makefile into 
      ../CubicIdeal.
   Compile libcubic via make in ../CubicIdeal.
   Download infrastructure_ideal.cc, infrastructure_ideal.h and makefile 
      into ../Infideal.
   Compile libinfideal via make in ../Infideal.
   Download hbounds.cc, hbounds.h and makefile 
      into ../Hbounds.
   Compile libhbounds via make in ../Hbounds.
   Download classnumber.cc, classnumber.h and makecn 
      into ../PointCount.
   Compile classnumber via ./makecn in ../PointCount.
   
   Make sure that the paths to the NTL libraries are
   specified correctly in the makefiles.
   
   Known Problems: 
   ---------------
   For unit rank 2 - Need to incorporate Kangaroo methods for extracting
                     R from h for larger examples. (Phase 4)
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   4. Incorporate Kangaroo methods into Phases 3 and 4 - Check Fontein's
      paper for improving Phase 4.
   5. Organize this file better and possibly even split it up into several
      files that are more easily managed.

   Author: 
   Eric Landquist
   
**classnumber.cc**

   This file contains functions to compute the divisor class number, h, 
   of any purely cubic function field. It also contains routines to       
   determine the regulator, R, of a cubic function field of positive unit 
   rank via computing in the infrastructure. Included are routines for
   applying the Baby Step-Giant Step algorithm of Shanks and also the 
   Kangaroo method of Pollard.
   For the functions, we will use results and algorithms from:

   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)
   
   
   How to create the executable:
   -----------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download cubic_ideal.cpp, cubic_ideal.h and makefile into 
      ../CubicIdeal.
   Compile libcubic via make in ../CubicIdeal.
   Download infrastructure_ideal.cc, infrastructure_ideal.h and makefile 
      into ../Infideal.
   Compile libinfideal via make in ../Infideal.
   Download hbounds.cc, hbounds.h and makefile 
      into ../Hbounds.
   Compile libhbounds via make in ../Hbounds.
   Download classnumber.cc, classnumber.h and makecn 
      into ../PointCount.
   Compile classnumber via ./makecn in ../PointCount.
   
   Make sure that the paths to the NTL libraries are
   specified correctly in the makefiles.
   
   Known Problems: 
   ---------------
   For unit rank 2 - Test Kangaroo method in Phases 3 and 4.
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   4. Organize this file better and possibly even split it up into several
      files that are more easily managed.

   Author: 
   Eric Landquist
