function [ xTF ] = DiffusionEquation( hho, source  )
%DIFFUSIONEQUATION Solves the diffusion equation with zero Dirichlet boundary
% conditions:
%
%         { -div( grad(u)) = f     in Omega
%         { u = 0     on Boundary(Omega)
%

ncell_dofs = hho.mesh.ncells * hho.elements{1}.ncell_dofs;
nedge_dofs = hho.mesh.nedges * hho.elements{1}.nedge_dofs;
ntotal_dofs = ncell_dofs + nedge_dofs;

A = sparse(nedge_dofs, nedge_dofs);
A_sc = sparse(ncell_dofs, nedge_dofs);

B = zeros(nedge_dofs, 1);
B_sc = zeros(ncell_dofs, 1);

for cell_i = 1:hho.mesh.ncells
    
    element = hho.elements{cell_i};
    dirichlet_edges = find(hho.mesh.cell_neighbors{cell_i}==0);
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
    % Local contribution into global matrix
    for i = 1:element.nedges
        i_global = hho.mesh.cell_edges{cell_i}(i);
        for ik = 1:element.nedge_dofs
            i_local = (i-1) * element.nedge_dofs + ik;
            ipos = (i_global - 1) * element.nedge_dofs + ik;
            if all(dirichlet_edges~=i) % if non-dirichlet edge
                for j = 1:element.nedges
                    j_global = hho.mesh.cell_edges{cell_i}(j);
                    for jk = 1:element.nedge_dofs
                        j_local = (j-1) * element.nedge_dofs + jk;
                        jpos = (j_global - 1) * element.nedge_dofs + jk;
                        A(ipos,jpos) = A(ipos,jpos) + AFF(i_local,j_local);
                    end
                end
                B(ipos) = B(ipos) + bF(i_local);
            else % dirichlet edge, impose dirichlet BC here
                A(ipos,ipos)=1;   %since homogeneous dirichlet BC, we know that u=0 on all boundary edges
            end
            
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
    
end

% Solve for the edge degrees of freedom
xF = A \ B;

% Recover cell degrees of freedom from static condensation
xTF = zeros(ntotal_dofs,1);
xTF(1:ncell_dofs) = B_sc - A_sc * xF(1:nedge_dofs);
xTF(1+ncell_dofs:ntotal_dofs) = xF;

end

