This branch contains project files for the poroelastic system (both coupled and uncoupled to a ODE system at the boundary). The schemes utilize finite difference method for obtaining solutions and are written in python. 

FiniteDiffBiot1-D.py is a finite difference solver for the 1-D Biot System with Dirichlet Boundary Conditions written in python under assumptions of neglible inertia, fully saturated mixtures and small deformations of the domain. 

FiniteDiffCoupled1-D.py is a finite difference solver for the 1-D coupled system between a scaler ODE and the poroelastic Biot system connected at the right most endpoint of an interval (0,L). The time discretization utilizes a backwards euler method. 
