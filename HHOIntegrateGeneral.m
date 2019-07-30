function val = HHOIntegrateGeneral( hho, f )
%HHOIntegrateGeneral Integrate the given function f using the HHO
% quadrature points. Accurate for polynomials up to order ~3K. The function
% f must take three parameters (T, x, y) where T is the cell id number
    val = 0;
    for cell_i = 1:hho.mesh.ncells
        element = hho.elements{cell_i};
        for edge_i = 1:element.nedges
            for iqn = 1:element.ncell_qnodes
                xQN = element.cell_quad_rules{edge_i}(iqn,1);
                yQN = element.cell_quad_rules{edge_i}(iqn,2);
                wQN = element.cell_quad_rules{edge_i}(iqn,3);
                val = val + wQN * f(cell_i, xQN, yQN);
            end
        end
    end
end

