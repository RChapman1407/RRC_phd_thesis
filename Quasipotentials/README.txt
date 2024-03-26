Here, I will briefly describe the changes which have been made to the code "olim4vad.c" to run for our system, 
as well as changes which must be made to run different calibrations and mesh refinements. 
Line numbers given correspond to those in "AMOC_OLIM4_decades_N_1024_K_22_FAMOUSB2xCO2_on_TOL1e-12.c"

Line 19, We set to the tolerance to be TOL  = 1.0 x 10^{-12}
Line 20 and 21, define the mesh spacing NX, NY
Line 22, define the step size K
Line 27, include the parameter file for the calibration, "FAMOUSB2xCO2_on.h"
Lines 158-164, Names of output files can be changed. Vector fields are also being output.
Lines 168-201,  include the equations for the system being calculated
Lines 205-216, include the noise matrix, Sigma. Here, we set the noise to be the identity matrix, unless otherwise stated.
Line 253, include the Jacobian evaluated at the starting point
In the parameter file, Lines 39-42, define the boundaries of the calculation
In the parameter file, Lines 46 and 47, define the starting position for the calculation
In the parameter file, Lines 51 and 52, define the end point of a minimum action path for the backwards shooting method.
