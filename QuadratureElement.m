% ----------------------------------------------------------------------
%           Quadrature information for a cell of a 2D mesh
% ----------------------------------------------------------------------
%
% Author: Daniel Anderson
%
% The code for computing the Dunavant quadrature rules was downloaded
% from https://people.sc.fsu.edu/~jburkardt/m_src/triangle_dunavant_rule/triangle_dunavant_rule.html
% ----------------------------------------------------------------------
classdef QuadratureElement
    %QUADRATUREELEMENT Data structure containing quadrature rules on a cell
    %   Stores the quadrature rules for a cell and its adjacent edges, as
    %   well as evaluations of the basis functions at each quadrature node.
    
    properties
        cell_id                 % ID number of the corresponding cell
        nedges                  % The number of edges of this element
        cell_doe                % The degree-of-exactness of the cells
        edge_doe                % The degree-of-exactness on the edges
        cell_quad_rules         % Quadrature rules on the cells
        edge_quad_rules         % Quadrature rules on the edges
        ncell_qnodes            % Number of quadrature nodes for cells
        nedge_qnodes            % Number of quadrature nodes for edges
        cell_values_in_cell     % Values of the cell basis functions
        cell_values_on_edge     % Values of the cell basis functions
        cell_gradients_in_cell  % Gradients of the basis inside the cell
        cell_gradients_on_edge  % Gradients of the basis on the edges
        edge_values             % Values of the edge basis functions
        nhighorder_dofs         % Number of K+1 cell degrees of freedom
        ncell_dofs              % Number of cell degrees of freedom
        nedge_dofs              % Number of edge degrees of freedom
        ngrad_dofs              % Dimension of the gradient space
        ntotal_edge_dofs        % Total number of edge degrees of freedom
        ntotal_dofs             % Total number of degrees of freedom
    end
    
    methods
        function obj = QuadratureElement(cell_id,vertices,center,K,cell_doe,edge_doe,cell_basis,edge_bases)
            %QUADRATUREELEMENT Construct a quadrature element on the given cell
            %   cell_id     The index of the cell
            %   vertices    The vertices of the cell as an n x 2 matrix
            %   center      The cell center point as a 1 x 2 matrix
            %   K           The degree of the cell polynomials
            %   cell_doe    The degree of exactness for cell quadrature
            %   edge_doe    The degree of exactness for edge quadrature
            %   cell_basis  The basis functions for polynomials on the cell. These should be given as an n x 1 matrix of CellBasisFunction objects
            %   edge_basis  The basis functions for polynomials on the edges. These should be given as an n x 1 matrix of EdgeBasisFunction objects
            
            obj.cell_id = cell_id;
            obj.nedges = length(vertices)-1;
            obj.cell_doe = cell_doe;
            obj.edge_doe = edge_doe;
            
            obj.nhighorder_dofs = (K+2)*(K+3)/2;
            obj.ncell_dofs = (K+1)*(K+2)/2;
            obj.ngrad_dofs = (K+2)*(K+3)/2 - 1;
            obj.nedge_dofs = length(edge_bases{1});
            obj.ntotal_edge_dofs = obj.nedges * obj.nedge_dofs;
            obj.ntotal_dofs = obj.ncell_dofs + obj.ntotal_edge_dofs;
            
            % Compute quadrature rules for the cell domains
            obj.cell_quad_rules = cell(obj.nedges,1);
            for edge_i = 1:obj.nedges
                triangle = [vertices(edge_i,:); vertices(edge_i+1,:); center];
                obj.cell_quad_rules{edge_i} = obj.getDunavantRule(triangle,cell_doe);
            end
            
            % Compute quadrature rules for the edge domains
            obj.edge_quad_rules = cell(obj.nedges,1);
            for edge_i = 1:obj.nedges
                interval = [vertices(edge_i,:); vertices(edge_i+1,:)];
                obj.edge_quad_rules{edge_i} = obj.getGaussianRule(interval,edge_doe);
            end
            
            % Count quadrature nodes
            obj.ncell_qnodes = size(obj.cell_quad_rules{1},1);
            obj.nedge_qnodes = size(obj.edge_quad_rules{1},1);
            
            % Compute quadrature for cell basis functions on the cell
            obj.cell_values_in_cell = zeros(obj.nedges, length(cell_basis), obj.ncell_qnodes);
            for edge_i = 1:obj.nedges
                for b = 1:length(cell_basis)
                    for iqn = 1:obj.ncell_qnodes
                        pt = obj.cell_quad_rules{edge_i}(iqn, 1:2);
                        obj.cell_values_in_cell(edge_i,b,iqn) = cell_basis{b}.eval(pt(1), pt(2));
                    end
                end
            end
            
            % Compute quadrature for cell basis functions on the edges
            obj.cell_values_on_edge = zeros(obj.nedges, length(cell_basis), obj.nedge_qnodes);
            for edge_i = 1:obj.nedges
                for b = 1:length(cell_basis)
                    for iqn = 1:obj.nedge_qnodes
                        pt = obj.edge_quad_rules{edge_i}(iqn, 1:2);
                        obj.cell_values_on_edge(edge_i,b,iqn) = cell_basis{b}.eval(pt(1), pt(2));
                    end
                end
            end
            
            % Compute quadrature for gradients in the cells
            obj.cell_gradients_in_cell = zeros(obj.nedges, length(cell_basis), obj.ncell_qnodes, 2);
            for edge_i = 1:obj.nedges
                for b = 1:length(cell_basis)
                    for iqn = 1:obj.ncell_qnodes
                        pt = obj.cell_quad_rules{edge_i}(iqn, 1:2);
                        obj.cell_gradients_in_cell(edge_i,b,iqn,:) = cell_basis{b}.eval_grad(pt(1), pt(2));
                    end
                end
            end
            
            % Compute quadrature for gradients on the edges
            obj.cell_gradients_on_edge = zeros(obj.nedges, length(cell_basis), obj.nedge_qnodes, 2);
            for edge_i = 1:obj.nedges
                for b = 1:length(cell_basis)
                    for iqn = 1:obj.nedge_qnodes
                        pt = obj.edge_quad_rules{edge_i}(iqn, 1:2);
                        obj.cell_gradients_on_edge(edge_i,b,iqn,:) = cell_basis{b}.eval_grad(pt(1), pt(2));
                    end
                end
            end
            
            % Compute quadrature for edge basis functions
            obj.edge_values = zeros(obj.nedges, length(edge_bases{1}), obj.nedge_qnodes);
            for edge_i = 1:obj.nedges
                for b = 1:length(edge_bases{edge_i})
                    for iqn = 1:obj.nedge_qnodes
                        pt = obj.edge_quad_rules{edge_i}(iqn, 1:2);
                        obj.edge_values(edge_i,b,iqn) = edge_bases{edge_i}{b}.eval(pt(1), pt(2));
                    end
                end
            end
        end
        
        function qr = getDunavantRule(~,triangle,doe)
            %GetDunavantRule Compute a Dunavant quadrature rule for
            % quadrature on the triangular subelements of a cell
            
            % The dunavant rule gives quadrature for the unit triangle,
            % so we scale and transform the given quadrature rule for
            % our triangle in question
            measure = triangle_area(triangle');
            [xy, w] = dunavant_rule(doe);
            qr = zeros(length(w),3);
            
            c1 = triangle(1,:) - triangle(3,:);
            c2 = triangle(2,:) - triangle(3,:);
            M = [c1' c2'];
            
            % Scale and translate the quadrature rule
            for i = 1:length(w)
                qr(i,1:2) = (M * xy(:,i) + triangle(3,:)')';
                qr(i,3) = w(i) * measure;
            end
        end   
        
        function qr = getGaussianRule(~,interval,doe)
            %GetGaussianRule Compute a Gaussian quadrature rule for
            % quadrature an edge (an interval embedded in 2D)
            
            % Compute number of nodes for the desired degree of exactness
            if mod(doe, 2) == 0; doe = doe + 1; end
            num_nodes = (doe + 1) / 2;
            
            % Compute Gaussian quadrature nodes
            if num_nodes == 1
                weights = [1];
                nodes = [0.5];
            else
                M = zeros(num_nodes, num_nodes);
                for i = 2:num_nodes
                    M(i, i-1) = sqrt(1/(4 - 1/(i-1)^2));
                    M(i-1, i) = M(i, i-1);
                end
                [V,D] = eig(M);
                weights = (V(1,:).^2)';
                nodes = (diag(D)+1)/2;
            end
            
            % Scale the nodes since the quadrature rule is over [0,1]
            measure = norm(interval(2,:) - interval(1,:));
            qr = zeros(num_nodes,3);

            for i=1:num_nodes
               qr(i,1:2) = interval(1,:) + nodes(i) * (interval(2,:) - interval(1,:));
               qr(i,3) = weights(i) * measure;
            end     
        end
        
    end
end

