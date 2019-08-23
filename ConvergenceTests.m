% Check the convergence of the HHO diffusion scheme

% If you use this code or part of it in a scientific publication, 
% please mention the following book as a reference for the underlying principles
% of HHO schemes:
%
% The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
% D. A. Di Pietro and J. Droniou. 2019, 516p. 
% url: https://hal.archives-ouvertes.fr/hal-02151813.


% Change this to the directory containing the meshes
mesh_directory = 'matlab_meshes/';

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
    {'mesh2_1', 'mesh2_2', 'mesh2_3'}%, 'mesh2_4', 'mesh2_5'}
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
% The exact solution, its gradient, 
% and the source term f = -div(\nabla u)

    u_exact_D =  @(x,y) sin(pi*x)*sin(pi*y);
    gradu_exact = @(x,y) [pi*cos(pi*x)*sin(pi*y); pi*sin(pi*x)*cos(pi*y)];
    source = @(x,y) 2*pi^2*sin(pi*x)*sin(pi*y);
    norm_gradu = @(T,x,y) norm(gradu_exact(x,y))^2;

% -------------------------------------------------------------------
% Run tests 

clf;

h_size = cell(length(mesh_names), MAX_K);
error_grad = cell(length(mesh_names), MAX_K);
error_D = cell(length(mesh_names), MAX_K);

for mesh_i = 1:length(mesh_names)
    fprintf('*** %s ***\n', mesh_names{mesh_i});
    
    % Compute relative errors and mesh sizes
    for K = 0:MAX_K    
        fprintf('K = %d::\n', K);
        h_size{mesh_i,K+1} = zeros(length(meshes{mesh_i}), 1);
        error_grad{mesh_i,K+1} = zeros(length(meshes{mesh_i}), 1);
        error_D{mesh_i,K+1} = zeros(length(meshes{mesh_i}), 1);
        for mesh_j = 1:length(meshes{mesh_i})
            hho = HHO(strcat(mesh_directory, meshes{mesh_i}{mesh_j}, '.mat'), K);
            h_size{mesh_i,K+1}(mesh_j) = max(hho.mesh.h_size);
                u_D = DiffusionEquation(hho,source);
                [grad_pT]  = coeffGradP( hho, u_D);
                error_D{mesh_i,K+1}(mesh_j) = HHORelError(hho, u_D, u_exact_D);
                fprintf('\t%s, h = %f, error_u = %f\n', meshes{mesh_i}{mesh_j}, h_size{mesh_i,K+1}(mesh_j), error_D{mesh_i,K+1}(mesh_j));
                diff = @(T,x,y) norm(GradP(hho,grad_pT,T,x,y)-gradu_exact(x,y))^2;
            error_grad{mesh_i,K+1}(mesh_j) = sqrt(HHOIntegrateGeneral(hho, diff))/sqrt(HHOIntegrateGeneral(hho,norm_gradu));
            fprintf('\t%s, h = %f, error_gradu = %f\n', meshes{mesh_i}{mesh_j}, h_size{mesh_i,K+1}(mesh_j), error_grad{mesh_i,K+1}(mesh_j));
        end
    end
    
    % Plot convergence figure for this mesh family
    
    figure(mesh_i);
    axes('XScale', 'log', 'YScale', 'log')
    xlim auto;
    ylim auto;
    hold on;
    title(strcat(mesh_names{mesh_i}, ', Dirichlet BC'));
    xlabel('h');
    ylabel('L^2 Error, u');
    
    for K=0:MAX_K
        loglog(h_size{mesh_i,K+1}, error_D{mesh_i,K+1}, 'DisplayName', strcat('K = ', num2str(K)));
    end
    
    figure(mesh_i+1);
    axes('XScale', 'log', 'YScale', 'log')
    xlim auto;
    ylim auto;
    hold on;
    title(strcat(mesh_names{mesh_i}, ', Dirichlet BC'));
    xlabel('h');
    ylabel('L^2 Error, gradu');
    
    for K=0:MAX_K
        loglog(h_size{mesh_i,K+1}, error_grad{mesh_i,K+1}, 'DisplayName', strcat('K = ', num2str(K)));
    end

end
