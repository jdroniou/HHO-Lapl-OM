% Example usage of the HHO diffusion equation solver
% Change this to the directory containing the meshes
mesh_directory = 'matlab_meshes/';

% -------------------------------------------------------------------
% Setup problem parameters
% Exact solution -
u_exact = @(x,y) sin(pi*x)*sin(pi*y);
% Source term - takes parameters T,x,y where T is the cell number
source = @(x,y) 2*pi^2*sin(pi*x)*sin(pi*y);
% -------------------------------------------------------------------
% Construct the HHO data structure

% mesh = 'mesh1_3.mat';     % triangular mesh
% mesh = 'hexa1_2.mat';     % hexa mesh
mesh = 'mesh2_3.mat';       % Cartesian mesh
% mesh = 'mesh4_1.mat';     % Kershaw mesh
K = 1;                      % The polynomial degree on the cells and edges

hho = HHO(strcat(mesh_directory, mesh), K);

% -------------------------------------------------------------------
% Solve the problem

u = DiffusionEquation(hho, source);

% -------------------------------------------------------------------
% Plot the solution
error = HHORelError(hho, u, u_exact);
fprintf('\nError in L2 norm: %d\n',error)

[ucell,uedge] = HHO_Cell_Edge_Ave(hho,u); %compute the average value of u on the cell and on the edges(for plotting)
uVal = [ucell;uedge];
ncell = hho.mesh.ncells;
nedge = hho.mesh.nedges;
nvert = hho.mesh.nverts;
cell_v = hho.mesh.cell_vertices;
cell_n = hho.mesh.cell_neighbors;
cell_e = hho.mesh.cell_edges;
vertex = hho.mesh.vertices;
write_solution_vtk(uVal,'solution_HHO',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

% -------------------------------------------------------------------
% Construct the gradient
grad_p = zeros(hho.mesh.ncells,2);
[grad_pT]  = coeffGradP( hho, u);
for cell_i=1:hho.mesh.ncells
    x = hho.mesh.cell_center(cell_i,1);
    y = hho.mesh.cell_center(cell_i,2);
    [grad_p(cell_i,:)]  = GradP( hho, grad_pT, cell_i, x, y); %compute the gradient at the cell centers
end
