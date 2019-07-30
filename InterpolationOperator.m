function [ a, b ] = InterpolationOperator( hho, cell_i, f )
%INTERPOLATIONOPERATOR The local interpolation operator

    element = hho.elements{cell_i};
    
    a = zeros(element.ntotal_dofs);
    b = zeros(element.ntotal_dofs, 1);

    % Compute quadrature over each edge of the cell
    for edge_i = 1:element.nedges
        % Compute volumetric terms
        for iqn = 1:element.ncell_qnodes
            xQN = element.cell_quad_rules{edge_i}(iqn,1);
            yQN = element.cell_quad_rules{edge_i}(iqn,2);
            wQN = element.cell_quad_rules{edge_i}(iqn,3);
            for i = 1:element.ncell_dofs
               phi_i = element.cell_values_in_cell(edge_i, i, iqn);
               for j = 1:element.ncell_dofs
                  phi_j = element.cell_values_in_cell(edge_i, j, iqn);
                  a(i, j) = a(i, j) + wQN * phi_i * phi_j;
               end
               b(i) = b(i) + wQN * phi_i * f(xQN, yQN);
            end
        end
        % Boundary terms
        for iqn = 1:element.nedge_qnodes
            xQN = element.edge_quad_rules{edge_i}(iqn,1);
            yQN = element.edge_quad_rules{edge_i}(iqn,2);
            wQN = element.edge_quad_rules{edge_i}(iqn,3);
            offset = element.ncell_dofs + (edge_i-1) * element.nedge_dofs;
            for i = 1:element.nedge_dofs
                phi_i = element.edge_values(edge_i, i, iqn);
                for j = 1:element.nedge_dofs
                     phi_j = element.edge_values(edge_i, j, iqn);
                     a(offset+i,offset+j) = a(offset+i,offset+j) + wQN * phi_i * phi_j;
                end
                b(offset+i) = b(offset+i) + wQN * phi_i * f(xQN, yQN);
            end
        end
    end
end

