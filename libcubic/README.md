**cubic_ideal.cpp** contains functions to perform ideal arithmetic such as 
   multiplication, division, squaring, inversion and reduction. All the 
   special cases have also been implemented with the addition of 
   functions to calculate the minimal element and the canonical basis. 
   For the functions, we will use results and algorithms from:

   [S01] R. Scheidler, "Ideal arithmetic and infrastructure in purely 
                        cubic function fields"
   [B04] M. Bauer, "The Arithmetic of Certain Cubic Function Fields"
   [B05] M. Bauer, Untitled Manuscript
   [L09] E. Landquist, "Infrastructure, Arithmetic, and Class Number
                        Computation in Cubic Function Fields of 
			Characteristic at Least 5" (PhD Thesis)

   There are also miscellaneous functions to calculate the gcd of three 
   elements and bubble sort which were implemented when required by the 
   specific functions. 
   
   How to create the library:
   --------------------------
   Download invariants.cc, invariants.h and makefile into ../Invariants.
   Compile libinvariants via make in ../Invariants.
   Download cubic_ideal.cpp, cubic_ideal.h and makefile into 
      ../CubicIdeal.
   Compile libcubic via make in ../CubicIdeal.
   
   Make sure that the paths to the NTL libraries are
   specified correctly in the makefile. 
   
   Known Problems: 
   ---------------
   None.
   
   Work to be completed:
   ---------------------
   1. Expand (via templates or other) to include zz_pX, zz_pEX, ZZ_pEX,  
      GF2X, and GF2EX. 
   2. Include support for characteristic 3 (p=3). See Jonathan Webster.
   3. Include explicit formulas where available: Flon and Oyono, "Fast 
      arithmetic on Jacobians of Picard curves" (deg(G) = 4, H=1)
   4. Include functionality for cubic number fields.
   5. Optimize where possible.

   Authors: 
   Mark Bauer
   Eric Nosal
   Eric Landquist
   Hameeduz Zaman, hameed_uz@hotmail.com, (403) 399-8582
