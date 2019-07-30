HHO scheme for the Poisson problem with homogeneous Dirichlet BC:

Main files: 

1.Example - runs an HHO scheme for the Poisson problem, computes and plots the solution u of an HHO scheme with order K, also measures the relative error;  	    
	    the gradient is also reconstructed, and its values on the cell centers are computed.
 
2.ConvergenceTests - runs a HHO scheme for the Poisson problem on a family of meshes; computes the relative errors in the solution and the reconstructed gradient 
		     for each mesh size and each K, in order to obtain the rate of convergence

Initialisation: (author: Daniel Anderson) 

1. HHO - creates an HHO data structure for a given degree K
2. QuadratureElement- sets up the quadrature rules on each cell and each edge
3. Mesh2D - loads the mesh structure
4. CellBasisFunction - basis functions in the cells
5. EdgeBasisFunction - basis functions on the edges
6. D2BasisFunction - 2 dimensional basis functions (to be used for the basis of the gradient)

Functions involved in setting up the Dunavant quadrature rule: (Author: John Burkardt)

1. dunavant_rule - returns the points and weights of a Dunavant rule
2. dunavant_suborder - the suborders for a Dunavant rule
3. dunavant_suborder_num -  number of suborders for a Dunavant rule
4. dunavant_subrule - a compressed Dunavant rule
5. i4_modp - returns the nonnegative remainder of I4 division
6. i4_wrap - forces an integer to lie between given limits by wrapping
7. triangle_area - returns area of a triangle

Setting up and solving the equation:

1. DiffusionOperator - local operator for the diffusion 
2. LoadOperator - operator for the source term
3. DiffusionEquation - sets up the global equation, and solves for the unknowns using static condensation  											      (Note: For nonhomogeneous Dirichlet BC, modify the else component of the "if" statement in lines 50-52) 

Reconstruction of the gradient:

1. coeffGradP - computes the coefficients for the reconstructed gradient
2. GradP - function that returns the value of the reconstructed gradient evaluated at a point (x_i,y_i) located in a cell cell_i.

Processing:

1. HHOInterpolate - projects a given function to the space of HHO functions of degree k (piecewise polynomials of degree k on cells and edges)
2. HHOIntegrateGeneral - integrates a function f over the domain omega by using the HHO quadrature points (accurate for polynomials up to degree 3K)
3. HHORelError - computes the relative error of a solution u obtained from an HHO scheme 
4. HHO_Cell_Edge_Ave - computes the average value of the solution on the cells and on the edges
5. write_solution_vtk - writes a vtk file for visualising the solution profile

The files in the folder Neumann BC enables the user to run test cases with Neumann BC (homogeneous or nonhomogeneous)