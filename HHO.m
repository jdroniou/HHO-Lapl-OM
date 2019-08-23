% If you use this code or part of it in a scientific publication, 
% please mention the following book as a reference for the underlying principles
% of HHO schemes:
%
% The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
% D. A. Di Pietro and J. Droniou. 2019, 516p. 
% url: https://hal.archives-ouvertes.fr/hal-02151813.


% ----------------------------------------------------------------------
%                 Core data structure for HHO schemes in 2D
% ----------------------------------------------------------------------
%
% Author: Daniel Anderson
%
% ----------------------------------------------------------------------
classdef HHO
    %HHO Data structure for HHO schemes
    %   This data structure contains all of the necessesities for
    %   implementing Hybrid High-Order schemes for generic 2D meshes. In
    %   particular, it curates the mesh, the basis functions and the
    %   quadrature rules for numerically integrating on cells and edges.
    
    properties
        mesh            % The mesh used by the scheme
        K               % The polynomial degree on the cells
        cell_bases      % The basis functions of the cells
        edge_bases      % The basis functions of the edges
        elements        % The quadrature elements for each mesh cell
    end
    
    methods
        function obj = HHO(mesh_filename, K)
            %HHO Construct an instance of this class
            %   Detailed explanation goes here
            obj.mesh = Mesh2D(mesh_filename);
            obj.K = K;
            
            % Degree of exactness for polynomial quadrature
            cell_doe = max(1, 3*K);
            edge_doe = max(1, 3*K);
            
            % Construct cell basis functions
            obj.cell_bases = cell(obj.mesh.ncells,1);
            for cell_i = 1:obj.mesh.ncells
                obj.cell_bases{cell_i} = obj.createCellBasis(cell_i);
            end
        
            % Construct edge basis functions
            obj.edge_bases = cell(obj.mesh.nedges,1);
            for edge_i = 1:obj.mesh.nedges
                obj.edge_bases{edge_i} = obj.createEdgeBasis(edge_i);
            end
            
            % Construct the quadrature elements
            obj.elements = cell(obj.mesh.ncells,1);
            for cell_i = 1:obj.mesh.ncells
                vertices = (obj.mesh.vertices(obj.mesh.cell_vertices{cell_i},:));
                center = obj.mesh.cell_center(cell_i,:);
                edge_bases_for_cell = obj.edge_bases(obj.mesh.cell_edges{cell_i});
                obj.elements{cell_i} = QuadratureElement(cell_i,vertices,center,K,cell_doe,edge_doe,obj.cell_bases{cell_i},edge_bases_for_cell);
            end
        end
        
        function offset = getBasisDegreeOffset(~,k)
            %getBasisDegreeOffset The offset into the basis array
            % at which elements of degree k begin
            offset = (k)*(k+1)/2;
        end
        
        function basis = createCellBasis(obj,cell_i)
            %createCellBasis Create a basis for polynomial functions on the
            % cell with the given index.
            
            % We create the cell basis with degree K+1 for the high-order
            % reconstruction operators.
            basis = cell((obj.K+2)*(obj.K+3)/2, 1);
            n = 1;
            for i = 0:obj.K+1
                for j = 0:i
                    xT = obj.mesh.cell_center(cell_i, :);
                    hT = obj.mesh.h_size(cell_i);
                    basis{n} = CellBasisFunction(xT, hT, j, i-j);
                    n = n + 1;
                end
            end
        end
        
        function basis = createEdgeBasis(obj,edge_i)
            %createEdgeBasis Create a basis for polynomial functions on the
            % cell edge with the given index
            basis = cell(obj.K + 1, 1);
            n = 1;
            for i = 0:obj.K
                center = obj.mesh.edge_center(edge_i,:);
                vertex_ids = obj.mesh.edge_vertices(edge_i, :);
                vertex_a = obj.mesh.vertices(vertex_ids(1),:);
                len = obj.mesh.edge_length(edge_i);
                basis{n} = EdgeBasisFunction(center,vertex_a,len,i);
                n = n + 1;
            end
        end

    end
end

