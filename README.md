# Cubic-Function-Fields

**distribution** contains files to compute the divisor class number and regulator of purely cubic function fields of characteristic at least 5, along with various statistics about the distribution of the class number. The class numbers are computed using the baby-step giant-step algorithm and the kangaroo algorithm, both of which I adapted to these fields.

**hbounds** contains files to make the library libhbounds, which contains utilities for computing an estimate, E, of the divisor class number of the cubic function field.

**infideal** contains files to make the library libinfideal, which contains functions to perform infrastructure arithmetic such as baby steps, giant steps, and inversion. Also included are routines to reduce ideals and bases, along with basis conversion functions and a method to extract the regulator from the divisor class number.

**invariants** contains files to make the library libinvariants, which contains global variables for use by any program computing in a cubic function field.

**kangaroo** contains client and server programs to compute the divisor class number of a large cubic function field using the kangaroo algorithm on a cluster or a computer with multiple processors.

**libcubic** contains files to make the library libcubic.a, which is a library to perform various arithmetic operations in the ideal class group and infrastructure of a cubic function field over a finite field of characteristic at least 5.

**landquist-thesis.pdf** is my Ph.D. dissertation, which contains the theory behind much of the code in this repository, along with the results of many of the calculations. This was submitted January 27, 2009.

**defense.pdf** is the slides of my Ph.D. defense, January 27, 2009.
