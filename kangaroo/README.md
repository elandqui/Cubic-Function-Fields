# Kangaroo algorithm

**kangaroo.cc** runs the individual tame and wild kangaroos for the parallelized kangaroo algorithm. This function makes the kangaroo jumps and checks for distinguished points, writing them to a file. The function terminates when the master program writes a file containing the solution. 

**rooking.cc** (King Kangaroo) is the master program for parallelized kangaroos. This program takes in information about a cubic function field, computes an approximation for its class number, an error bound, and information for the tame and wild kangaroo herds. It periodically opens the distinguished point file and checks for a match between the paths of a tame kangaroo and a wild kangaroo.
