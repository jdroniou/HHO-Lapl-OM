function  [grad_pT]  = coeffGradP( hho, p)
% grad_pT gives the coefficients for the polynomial that 
% reconstructs the gradient of an HHO solution to the Poisson equation
%
% hho      The HHO data structure
% p        The HHO solution to the Poisson equation
%
    grad_pT = cell(hho.mesh.ncells, 1);

    % Reconstruct the velocity in each cell
    for cell_i = 1:hho.mesh.ncells
        
        element = hho.elements{cell_i};
        
        % Compute the local gradient reconstruction operator
        [~, G_KT] = DiffusionOperator(hho, cell_i);
        
        % Compute local pressure in the cell
        pT = zeros(element.ntotal_dofs, 1);
        
        % Extract cell degrees of freedom
        for i = 1:element.ncell_dofs
           pT(i) = p((cell_i-1)*element.ncell_dofs + i);
        end
        
        % Extract edge degrees of freedom
        for j = 1:element.nedges
            j_global = hho.mesh.cell_edges{cell_i}(j);
            for jk = 1:element.nedge_dofs
                jpos = element.ncell_dofs*hho.mesh.ncells + (j_global - 1) * element.nedge_dofs + jk;
                j_local = element.ncell_dofs + (j-1) * element.nedge_dofs + jk;
                pT(j_local) = p(jpos);
            end
        end
        
        % Compute grad(p) in the cell
        grad_pT{cell_i} = G_KT * pT;
    end

end

