# libinfideal library

**infrastructure_ideal.cc**

   This file contains functions to perform infrastructure arithmetic such as 
   baby steps, giant steps, and inversion. Also included are routines to 
   reduce ideals and bases, along with basis conversion functions and a
   method to extract the regulator from a multiple. 
   For the functions, we will use results and algorithms from:

   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)
   [LSY] Y. Lee, R. Scheidler, and C. Yarrish, "Computation of the 
                        fundamental units and the regulator of a cyclic 
			cubic function field"
   [S01] R. Scheidler, "Ideal arithmetic and infrastructure in purely 
                        cubic function fields"
   
   How to create the library:
   --------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download cubic_ideal.cpp, cubic_ideal.h and makefile into 
      ../CubicIdeal.
   Compile libcubic via make in ../CubicIdeal.
   Download infrastructure_ideal.cc, infrastructure_ideal.h and makefile 
      into ../Infideal.
   Compile libinfideal via make in ../Infideal.

   Make sure that the paths to the NTL libraries are
   specified correctly in the makefile.
   
   Known Problems: 
   ---------------
   In below() for unit rank 2, there is an upper bound on y2 beyond which
   the function produces incorrect output (or perhaps none at all).
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, and ZZ_pEX.  
   2. Include functionality for cubic number fields.
   3. Optimize where possible.
   4. Fix below() to work more generally in unit rank 2.

   Authors: 
   Eric Landquist
   Christopher Yarrish 
