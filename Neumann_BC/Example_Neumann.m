% Example usage of the HHO diffusion equation solver
% Change this to the directory containing the meshes
mesh_directory = '\\ad.monash.edu\home\User027\hanzmarc\Documents\MATLAB\matlab_meshes\';

% -------------------------------------------------------------------
% Setup problem parameters
% -------------------------------------------------------------------
% Gives the exact solution, its gradient, 
% and the source term f = -div(\nabla u)
% Note: pure Neumann BC, so the solution computed is the one
% with 0 average, hence test cases here were chosen to be u such that u has
% zero average over the square (0,1)x(0,1).
global testcase;
testcase = 2;
if testcase == 1
    u_exact = @(x,y) cos(pi*x)*cos(pi*y);
    gradu_exact = @(x,y) [-pi*sin(pi*x)*cos(pi*y); -pi*cos(pi*x)*sin(pi*y)];
    source = @(x,y) 2*pi^2*cos(pi*x)*cos(pi*y);
elseif testcase == 2
    u_exact = @(x,y) sin(pi*x)*sin(pi*y)-4/pi^2;
    gradu_exact = @(x,y) [pi*cos(pi*x)*sin(pi*y); pi*sin(pi*x)*cos(pi*y)];
    source = @(x,y) 2*pi^2*sin(pi*x)*sin(pi*y);
elseif testcase == 3
    u_exact = @(x,y) cos(pi*x)*cos(pi*y)+sin(pi*x)*sin(pi*y)-4/pi^2;
    gradu_exact = @(x,y) [-pi*sin(pi*x)*cos(pi*y)+pi*cos(pi*x)*sin(pi*y); -pi*cos(pi*x)*sin(pi*y)+pi*sin(pi*x)*cos(pi*y)];
    source = @(x,y) 2*pi^2*cos(pi*x)*cos(pi*y)+2*pi^2*sin(pi*x)*sin(pi*y);
end

% -------------------------------------------------------------------
% Construct the HHO data structure

% mesh = 'mesh1_3';     % triangular mesh
% mesh = 'hexa1_2';     % hexa mesh
mesh = 'mesh2_3';       % Cartesian mesh
% mesh = 'mesh4_1';     % Kershaw mesh
K = 2;                      % The polynomial degree on the cells and edges

hho = HHO(strcat(mesh_directory, mesh), K);

% -------------------------------------------------------------------
% Solve the problem

[u] = DiffusionEquation_Neumann(hho, source, gradu_exact);

% -------------------------------------------------------------------
% Plot the solution
error = HHORelError(hho, u, u_exact);

[ucell,uedge] = HHO_Cell_Edge_Ave(hho,u); %compute the average value of u on the cell and on the edges(for plotting)
uVal = [ucell;uedge];
h_size = max(hho.mesh.h_size);
ncell = hho.mesh.ncells;
nedge = hho.mesh.nedges;
nvert = hho.mesh.nverts;
cell_v = hho.mesh.cell_vertices;
cell_n = hho.mesh.cell_neighbors;
cell_e = hho.mesh.cell_edges;
vertex = hho.mesh.vertices;
write_solution_vtk(uVal,'solution_HHO',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

% -------------------------------------------------------------------
% Construct the gradient and evaluate at cell centers
grad_p = zeros(hho.mesh.ncells,2);
[grad_pT]  = coeffGradP( hho, u);
for cell_i=1:hho.mesh.ncells
    x = hho.mesh.cell_center(cell_i,1);
    y = hho.mesh.cell_center(cell_i,2);
    [grad_p(cell_i,:)]  = GradP( hho, grad_pT, cell_i, x, y); %compute the gradient at the cell centers
end

% Measure the rel. error in the reconstructed gradient
diff = @(T,x,y) norm(GradP(hho,grad_pT,T,x,y)-gradu_exact(x,y))^2;
error_grad = sqrt(HHOIntegrateGeneral(hho, diff))/sqrt(HHOIntegrateGeneral(hho,norm_gradu));
            
fprintf('\t%s, h = %f, error_u = %f, error_gradu = %f \n', mesh, h_size , error, error_grad);
