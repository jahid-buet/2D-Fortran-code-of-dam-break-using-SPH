# 2D-Fortran-code-of-dam-break-using-SPH
This is a simple 2d code of dam break using Smoothed particles hydrodynamics(SPH) numerical method written in Fortran.
To run the simulation one has to run the executable file  "Dam_break.exe".
if you find any bug or problem in running the simulation u are free to change the source code situated in "Source files" folder.
To modify the source code one has to recompile the source files in "Source" folder to create new executable file(.exe).The source files were compiled with ifort compiler.


Definition of source files are summarized as follow:

ComDense.f90: Compute rate of change of desity of fluid particles using continuity equation.

ComGravity.f90: Compute acceleration due to gravity for all fluid particles in domain.

Compress.f90: Compute pressure between fluid particles by using equation of state.

ComPressGrad.f90:Compute accleration due to pressure gradient for all particles(boundary+fluid) in domain.

ComViscosity.f90: Compute acceleration due to vsicosity term for all fluid particles. Depending on  slip or no-slip boundary condition,viscosity for boundary particles are computed.

Dam_break.f90: Main source file from where all subroutines are called  and user interface for input parametrs are defined .Also time integration is implemented in this source file.

KernalGradientCorrection.f90: Compute inverse matrix related to formulation of kernal gradient coorection which is used in only in viscosity computation.

Time.f90:This module computes time step needs for simulation although in this  2d dam break simulation this routine is not used.

Xsph.f90:This module computes Xsph to regularize the particles movement in SPH.it is not used in this simulation.

geometry.f90:This module  compute the required geometry(3 wall and 1 water column) according to the dimension are given in user interface.

initial.f90:This module is called from main source file before time initegration step and initialize variables(velocity,acceleration and pressure ) for all particles in domain.

kernal.f90: This module computes kernel function needs for SPH simulation

kgf_matrix.f90:This module computes 3 by 3 inverse matrix for kernal gradient free SPH formulation  and used in Orginal viscosity.

part.f90: This module define derived data types of particles and its various component.

strain.f90:This module computes strain rate for fluid particles which is used in orginal viscosity formnulation(direct discretiztion of 1st and 2nd order derivatives using SPH function and then kgf matrix are used to improve the accuracy of viscosity formulation).

var.f90: this module is used to compue density and mass of all particlesin domain.

vector.f90:this model contains different derived data types which is used by part.f90 module.

Dam_break.exe: this is an executable file wihcih is created after compiling all source file in ifort compiler.


