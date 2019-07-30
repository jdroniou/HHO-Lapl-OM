function [ b ] = LoadOperator_NeumannBC(hho, cell_i, gradu_exact)
%LOADOPERATOR Compute the load operator (source term) for
% a diffusion equation for the HHO method

element = hho.elements{cell_i};
b = zeros(element.ntotal_edge_dofs, 1);

for edge_i = 1:element.nedges
    if hho.mesh.cell_neighbors{cell_i}(edge_i)==0 % detects Neumann edge
        for iqn = 1:element.nedge_qnodes
            % The quadrature point (xQN,yQN) and the weight wQN
            xQN = element.edge_quad_rules{edge_i}(iqn,1);
            yQN = element.edge_quad_rules{edge_i}(iqn,2);
            wQN = element.edge_quad_rules{edge_i}(iqn,3);
            f_xQN = gradu_exact(xQN, yQN);
            for i = 1:element.nedge_dofs
                phi_i = element.edge_values(edge_i, i, iqn);
                b((edge_i-1)*element.nedge_dofs+i) = b((edge_i-1)*element.nedge_dofs+i) + ...
                    wQN * phi_i * hho.mesh.normal{cell_i}(edge_i,:)*(f_xQN);
            end
        end
    end
end
end