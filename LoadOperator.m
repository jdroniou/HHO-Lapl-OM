function [ b ] = LoadOperator(hho, cell_i, source)
%LOADOPERATOR Compute the load operator (source term) for
% a diffusion equation for the HHO method (in the cells)

    element = hho.elements{cell_i};
    b = zeros(element.ncell_dofs, 1);
    
    for edge_i = 1:element.nedges
        for iqn = 1:element.ncell_qnodes
            % The quadrature point (xQN,yQN) and the weight wQN
            xQN = element.cell_quad_rules{edge_i}(iqn,1);
            yQN = element.cell_quad_rules{edge_i}(iqn,2);
            wQN = element.cell_quad_rules{edge_i}(iqn,3);
            f_xQN = source( xQN, yQN);
            for i = 1:element.ncell_dofs
               phi_i = element.cell_values_in_cell(edge_i, i, iqn);
               b(i) = b(i) + wQN * phi_i * f_xQN;
            end
        end
    end
end

