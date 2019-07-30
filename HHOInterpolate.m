function u = HHOInterpolate( hho, f )
%HHOINTERPOLATE Produce an interpolated HHO function from a given function

    ncell_dofs = hho.mesh.ncells * hho.elements{1}.ncell_dofs;
    nedge_dofs = hho.mesh.nedges * hho.elements{1}.nedge_dofs;
    ntotal_dofs = ncell_dofs + nedge_dofs;

    u = zeros(ntotal_dofs, 1);
    
    for cell_i = 1:hho.mesh.ncells
        element = hho.elements{cell_i};
        % Compute local interpolation operator
        [a, b] = InterpolationOperator(hho, cell_i, f);
        UT = a \ b;
        % Cell terms
        u(1+(cell_i-1)*element.ncell_dofs:cell_i*element.ncell_dofs) = UT(1:element.ncell_dofs);
        % Edge terms
        for i = 1:element.nedges
            i_global = hho.mesh.cell_edges{cell_i}(i);
            cell_offset = hho.mesh.ncells * element.ncell_dofs;
            edge_offset = (i_global-1) * element.nedge_dofs;
            u(1+cell_offset+edge_offset:cell_offset+edge_offset+element.nedge_dofs) = ...
                UT(1 + element.ncell_dofs + (i-1)*element.nedge_dofs: element.ncell_dofs + i*element.nedge_dofs);
        end
    end     
end

