% ----------------------------------------------------------------------
%                   Data structure for 2D meshes
% ----------------------------------------------------------------------
%
% Author: Daniel Anderson
%
% ----------------------------------------------------------------------
classdef Mesh2D

    properties
        ncells              % The number of cells in the mesh
        nedges              % The number of edges in the mesh
        nbedges             % The number of boundary edges in the mesh
        nverts              % The number of vertices in the mesh
        ncell_edges         % The number of edges of each cell
        
        vertices            % The vertices of the cells of the mesh
        cell_vertices       % The indices of the vertices of each cell
        cell_edges          % The indices of the edges of each cell
        edge_vertices       % The indices of the verices of each edge
        
        cell_neighbors      % The neighbor(s) of each cell
        cell_center         % The midpoint of each cell
        edge_center         % The midpoint of each edge
        diameter            % The diameter of each cell
        edge_length         % The length of each edge
        area                % The area of each cell
        perimeter           % The perimeter of each cell
        h_size              % The ratio of each cell's area to perimeter
        normal              % The outer normals to the edges of the cells
    end
    
    methods
        
        function obj = Mesh2D(filename)
            %Mesh2D Load the mesh with the given filename.
            % The mesh file must be a MATLAB data file (.mat) in the
            % specified format given in matlab_meshes/README.txt
            
            raw_mesh = load(filename);
            
            % Read mesh size data
            obj.ncells = raw_mesh.ncell;
            obj.nedges = raw_mesh.nedge;
            obj.nverts = raw_mesh.nvert;
            obj.ncell_edges = zeros(obj.ncells,1);
            
            % Count number of boundary edges
            obj.nbedges = 0;
            for cell_i = 1:obj.ncells
               obj.nbedges = obj.nbedges + nnz(~raw_mesh.cell_n{cell_i});
            end
            
            % Read vertex and edge data
            obj.vertices = raw_mesh.vertex;
            obj.cell_vertices = raw_mesh.cell_v;
            obj.cell_edges = raw_mesh.cell_e;
            
            % Read precomputed statistics
            obj.cell_center = raw_mesh.center;
            obj.diameter = raw_mesh.diam;
            obj.area = raw_mesh.area;
            obj.cell_neighbors = raw_mesh.cell_n;
            
            % Populate sizes, edge vertices, centers, length and normals
            obj.edge_vertices = zeros(obj.nedges,2);
            obj.edge_center = zeros(obj.nedges,2);
            obj.edge_length = zeros(obj.nedges,1);
            obj.normal = cell(obj.ncells, 1);
            obj.perimeter = zeros(obj.ncells,1);
            obj.h_size = zeros(obj.ncells,1);
            
            for cell_i = 1:obj.ncells
                obj.ncell_edges(cell_i) = length(obj.cell_edges{cell_i});
                obj.normal{cell_i} = zeros(obj.ncell_edges(cell_i), 2);
                for edge_j = 1:length(obj.cell_edges{cell_i})
                    edge_id = obj.cell_edges{cell_i}(edge_j);
                    v = [obj.cell_vertices{cell_i}(edge_j), obj.cell_vertices{cell_i}(edge_j+1)];
                    vertex_a = obj.vertices(v(1),:);
                    vertex_b = obj.vertices(v(2),:);
                    obj.edge_vertices(edge_id,:) = v;
                    obj.edge_center(edge_id,:) = (vertex_a + vertex_b) / 2;
                    obj.edge_length(edge_id) = norm(vertex_b - vertex_a);
                    obj.normal{cell_i}(edge_j, :) = ([0 1; -1 0] * (vertex_b - vertex_a)' / norm(vertex_b - vertex_a))';
                    obj.perimeter(cell_i) = obj.perimeter(cell_i) + obj.edge_length(edge_id);
                end
                obj.h_size(cell_i) = obj.area(cell_i) / obj.perimeter(cell_i);
            end
        end
        
        % Determine which cell a given point is located in. Takes O(Nedges)
        % time per query, so use sparingly. If this is needed frequently
        % and speed is a concern, you should implement the proper algorithm
        % for point location, which takes O(log(Nedges)) per query, but
        % would probably take a year to implement in MATLAB.
        %
        % NOTE: If a point is contained in/on many cells (ie. it lies on an
        % edge or vertex), then one of them will be returned arbitrarily
        function cell_i = point_location(obj, x, y)
            for i = 1:obj.ncells
                if obj.point_in_cell(i, x, y)
                    cell_i = i;
                end
            end
        end
        
        % Determine whether a given point (x,y) lies in or on the cell
        % cell_i.
        function [in, on] = point_in_cell(obj, cell_i, x, y)
            xv = obj.vertices(obj.cell_vertices{cell_i}(1:end-1), 1);
            yv = obj.vertices(obj.cell_vertices{cell_i}(1:end-1), 2);
            [in, on] = inpolygon(x, y, xv, yv);
        end
    end
end

