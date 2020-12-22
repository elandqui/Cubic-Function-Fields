**everything.cc** Finds an approximation, E, of the divisor class number h and a bound L such that |h-E|<L^2 and gives statistics about the distribution.

**latticedata.cc** 

   This file contains functions to compute the divisor class number, h, 
   and regulator, R, of several purely cubic function fields of unit rank 2.
   We will also gather information about the fundamental units and the structure
   of the infrastructure. Included are routines for applying the Baby Step-Giant 
   Step algorithm of Shanks and also the Kangaroo method of Pollard.
   For the functions, we will use results and algorithms from:

   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)
   
   
   How to create the executable:
   --------------------------
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
   Compile latticedata via ./makeld in ../PointCount.
   
   Make sure that the paths to the NTL libraries are
   specified correctly in the makefiles.
   
   Known Problems: 
   ---------------
   Needs to be written.
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   4. Organize this file better and possibly even split it up into several
      files that are more easily managed.

   Author: 
   Eric Landquist
