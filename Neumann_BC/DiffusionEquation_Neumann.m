function [ xTF] = DiffusionEquation_Neumann( hho, source, gradu_exact  )
%DIFFUSIONEQUATION Solves the diffusion equation with Neumann boundary
% conditions:
%
%         { -div(K grad(u)) = f     in Omega
%         { K grad(u) * nTF = g     on Boundary(Omega)
%         { average(u)      = 0
%
% The solution normalised such that its average is zero.

    ncell_dofs = hho.mesh.ncells * hho.elements{1}.ncell_dofs;
    nedge_dofs = hho.mesh.nedges * hho.elements{1}.nedge_dofs;
    ntotal_dofs = ncell_dofs + nedge_dofs;
    
    nlocal_cell_dofs = hho.elements{1}.ncell_dofs;
    nlocal_edge_dofs = hho.elements{1}.nedge_dofs;

    A = sparse(nedge_dofs, nedge_dofs);
    A_sc = sparse(ncell_dofs, nedge_dofs);
    
    B = zeros(nedge_dofs, 1);
    B_sc = zeros(ncell_dofs, 1);
    
    cell_quadrature = zeros(ncell_dofs, 1);
    total_measure = 0;
    
    for cell_i = 1:hho.mesh.ncells
        mT = hho.mesh.area(cell_i);
        total_measure = total_measure + mT;
        
        element = hho.elements{cell_i};
        
        a = DiffusionOperator(hho, cell_i);
        b = LoadOperator(hho, cell_i, source);
        
        % Perform static condensation
        ATT = a(1:element.ncell_dofs, 1:element.ncell_dofs);
        ATF = a(1:element.ncell_dofs, element.ncell_dofs+1:element.ntotal_dofs);
        AFF = a(element.ncell_dofs+1:element.ntotal_dofs, element.ncell_dofs+1:element.ntotal_dofs);
        invATT_ATF = (ATT \ ATF);
        invATT_b = ATT \ b;
        AFF = AFF - ATF' * invATT_ATF;
        bF = -ATF' * invATT_b;
        bF = bF + LoadOperator_NeumannBC(hho, cell_i, gradu_exact);
        % Local contribution into global matrix
        for i = 1:element.nedges
            i_global = hho.mesh.cell_edges{cell_i}(i);
            for ik = 1:element.nedge_dofs
                i_local = (i-1) * element.nedge_dofs + ik;
                ipos = (i_global - 1) * element.nedge_dofs + ik;
                for j = 1:element.nedges
                    j_global = hho.mesh.cell_edges{cell_i}(j);
                    for jk = 1:element.nedge_dofs
                        j_local = (j-1) * element.nedge_dofs + jk;
                        jpos = (j_global - 1) * element.nedge_dofs + jk;
                        A(ipos,jpos) = A(ipos,jpos) + AFF(i_local,j_local);
                    end
                end
                B(ipos) = B(ipos) + bF(i_local);
            end
        end
        
        % Build static condensation operator
        B_sc(1+(cell_i-1)*element.ncell_dofs:cell_i*element.ncell_dofs) = invATT_b;
        for i = 1:element.ncell_dofs
            for j = 1:element.nedges
                j_global = hho.mesh.cell_edges{cell_i}(j);
                for jk = 1:element.nedge_dofs
                    ipos = (cell_i-1) * element.ncell_dofs + i;
                    jpos = (j_global - 1) * element.nedge_dofs + jk;
                    j_local = (j-1) * element.nedge_dofs + jk;
                    A_sc(ipos,jpos) = A_sc(ipos,jpos) + invATT_ATF(i, j_local);
                end
            end
        end
        
        % Assemble cell integrals
        for edge_i = 1:element.nedges
            for iqn = 1:element.ncell_qnodes
                wQN = element.cell_quad_rules{edge_i}(iqn,3);
                for i = 1:element.ncell_dofs
                   phi_i = element.cell_values_in_cell(edge_i, i, iqn);
                   ipos = (cell_i-1) * element.ncell_dofs + i;
                   cell_quadrature(ipos) = cell_quadrature(ipos) + wQN * phi_i;
                end
            end
        end
    end

    % Delete the first row to fix the unique solution
    A(1,:) = 0;
    A(1,1) = 1;
    B(1) = 0;

    % Solve for the edge degrees of freedom
    xF = A \ B;

    % Recover cell degrees of freedom from static condensation
    xTF = zeros(ntotal_dofs,1);
    xTF(1:ncell_dofs) = B_sc - A_sc * xF(1:nedge_dofs);
    xTF(1+ncell_dofs:ntotal_dofs) = xF;

    % Translate cells to zero average
    avg = cell_quadrature' * xTF(1:ncell_dofs) / total_measure;
%     avg = 0;
    for cell_i = 1:hho.mesh.ncells
        loc_i = 1 + (cell_i-1) * nlocal_cell_dofs;
        xTF(loc_i) = xTF(loc_i) - avg;
    end
    
    % Translate edges
    for edge_i = 1:hho.mesh.nedges
        loc_i = 1 + ncell_dofs + (edge_i-1) * nlocal_edge_dofs;
        xTF(loc_i) = xTF(loc_i) - avg;
    end
end

