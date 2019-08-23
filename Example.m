% Example usage of the HHO diffusion equation solver
% Change this to the directory containing the meshes
mesh_directory = 'matlab_meshes/';

% -------------------------------------------------------------------
% Setup problem parameters
% Exact solution -
u_exact = @(x,y) sin(pi*x)*sin(pi*y);
% gradient of exact solution
gradu_exact = @(x,y) [pi*cos(pi*x)*sin(pi*y); pi*sin(pi*x)*cos(pi*y)];
% Source term - 
source = @(x,y) 2*pi^2*sin(pi*x)*sin(pi*y);
% norm of the exact gradient
norm_gradu = @(T,x,y) norm(gradu_exact(x,y))^2;
% -------------------------------------------------------------------
% Construct the HHO data structure

% mesh = 'mesh1_3';     % triangular mesh
% mesh = 'hexa1_2';     % hexa mesh
mesh = 'mesh2_3';       % Cartesian mesh
% mesh = 'mesh4_1';     % Kershaw mesh
K = 1;                      % The polynomial degree on the cells and edges

hho = HHO(strcat(mesh_directory, mesh), K);

% -------------------------------------------------------------------
% Solve the problem

u = DiffusionEquation(hho, source);

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
% Construct the gradient and evaluate at cell centers, also compute rel.
% error for the reconstructed gradient
grad_p = zeros(hho.mesh.ncells,2);
[grad_pT]  = coeffGradP( hho, u);
for cell_i=1:hho.mesh.ncells
    x = hho.mesh.cell_center(cell_i,1);
    y = hho.mesh.cell_center(cell_i,2);
    [grad_p(cell_i,:)]  = GradP( hho, grad_pT, cell_i, x, y); %compute the gradient at the cell centers
end
diff = @(T,x,y) norm(GradP(hho,grad_pT,T,x,y)-gradu_exact(x,y))^2;
error_grad = sqrt(HHOIntegrateGeneral(hho, diff))/sqrt(HHOIntegrateGeneral(hho,norm_gradu));
            
fprintf('\t%s, h = %f, error_u = %f, error_gradu = %f \n', mesh, h_size , error, error_grad);
