% Check the convergence of the HHO diffusion scheme

% Change this to the directory containing the meshes
mesh_directory = '\\ad.monash.edu\home\User027\hanzmarc\Documents\MATLAB\matlab_meshes\';

% Maximum value of K to test (degree of the polynomial DOFs)
MAX_K = 1;

% Test meshes - Triangular, Cartesian, Kershaw, Hexagonal
mesh_choice = 2;
if mesh_choice == 1 % triangular
meshes = {
    {'mesh1_1', 'mesh1_2', 'mesh1_3', 'mesh1_4', 'mesh1_5'}
};
mesh_names = {'Triangular'};
elseif mesh_choice == 2 % Cartesian
meshes = {
    {'mesh2_1', 'mesh2_2', 'mesh2_3', 'mesh2_4', 'mesh2_5'}
};
mesh_names = {'Cartesian'};
elseif mesh_choice == 3 % Kershaw
meshes = {
    {'mesh4_1_1', 'mesh4_1_2', 'mesh4_1_3', 'mesh4_1_4', 'mesh4_1_5'}
};
mesh_names = {'Kershaw'};
elseif mesh_choice == 4 % hexagonal
  meshes = {
    {'hexa1_1', 'hexa1_2', 'hexa1_3', 'hexa1_4'}
};  
mesh_names = {'Hexagonal'};
else % do convergence test on all meshes
meshes = {
    {'mesh1_1', 'mesh1_2', 'mesh1_3', 'mesh1_4', 'mesh1_5'},...
    {'mesh2_1', 'mesh2_2', 'mesh2_3', 'mesh2_4', 'mesh2_5'},...
    {'mesh4_1_1', 'mesh4_1_2', 'mesh4_1_3', 'mesh4_1_4', 'mesh4_1_5'},...
    {'hexa1_1', 'hexa1_2', 'hexa1_3', 'hexa1_4'}
};
mesh_names = {'Triangular', 'Cartesian', 'Kershaw', 'Hexagonal'};
end
% -------------------------------------------------------------------
% Setup problem parameters

% -------------------------------------------------------------------
% Setup problem parameters
% -------------------------------------------------------------------
% Gives the exact solution, its gradient, 
% and the source term f = -div(\nabla u)
% Note: For pure Neumann BC, solution computed is the one
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
% Run tests (Neumann)

clf;

h_size = cell(length(mesh_names), MAX_K);
error = cell(length(mesh_names), MAX_K);

for mesh_i = 1:length(mesh_names)
    fprintf('*** %s ***\n', mesh_names{mesh_i});
    
    % Compute relative errors and mesh sizes
    for K = 0:MAX_K    
        fprintf('K = %d::\n', K);
        h_size{mesh_i,K+1} = zeros(length(meshes{mesh_i}), 1);
        error{mesh_i,K+1} = zeros(length(meshes{mesh_i}), 1);
        for mesh_j = 1:length(meshes{mesh_i})
            hho = HHO(strcat(mesh_directory, meshes{mesh_i}{mesh_j}, '.mat'), K);
            u = DiffusionEquation_Neumann(hho, source, gradu_exact);
            h_size{mesh_i,K+1}(mesh_j) = max(hho.mesh.h_size);
            error{mesh_i,K+1}(mesh_j) = HHORelError(hho, u, u_exact);
            fprintf('\t%s, h = %f, error = %f\n', meshes{mesh_i}{mesh_j}, h_size{mesh_i,K+1}(mesh_j), error{mesh_i,K+1}(mesh_j));
        end
    end
    
    % Plot convergence figure for this mesh family
    figure(mesh_i+1);
    axes('XScale', 'log', 'YScale', 'log')
    xlim auto;
    ylim auto;
    hold on;
    title(strcat(mesh_names{mesh_i}, ', Neumann BC'));
    xlabel('h');
    ylabel('L^2 Error');
    
    for K=0:MAX_K
        loglog(h_size{mesh_i,K+1}, error{mesh_i,K+1}, 'DisplayName', strcat('K = ', num2str(K)));
    end
    
   
end
