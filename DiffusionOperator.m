% ----------------------------------------------------------------------
%                        HHO Diffusion Operator
% ----------------------------------------------------------------------
%
% Author: Daniel Anderson
%
% ----------------------------------------------------------------------

function [A_TF, G_KT] = DiffusionOperator( hho, cell_i)
%DIFFUSIONOPERATOR Compute the local HHO diffusion operator a_T
%   hho         The HHO data structure
%   cell_i      The index of the cell T
% Returns the local diffusion operator matrix A_TF and the local
% gradient reconstruction operator matrix G_KT

    element = hho.elements{cell_i};
    
    grad_basis_offset = hho.getBasisDegreeOffset(1);

    % Initialise matrices
    MG = zeros(element.ngrad_dofs);
    BG = zeros(element.ngrad_dofs, element.ntotal_dofs);
    
    MTT = zeros(element.nhighorder_dofs);
    MFT = zeros(element.nedges, element.nedge_dofs, element.nhighorder_dofs);
    MFF = zeros(element.nedges, element.nedge_dofs, element.nedge_dofs);
    
    % Compute quadrature over each edge of the cell
    for edge_i = 1:element.nedges
        
        nTF = hho.mesh.normal{cell_i}(edge_i, :)';

        % Compute volumetric terms
        for iqn = 1:element.ncell_qnodes
            % The quadrature point (xQN,yQN) and the weight wQN
            xQN = element.cell_quad_rules{edge_i}(iqn,1);
            yQN = element.cell_quad_rules{edge_i}(iqn,2);
            wQN = element.cell_quad_rules{edge_i}(iqn,3);

            % Evaluate the diffusion tensor at this quadrature node
            kappa_iqn = eye(2,2);
            
            % Compute gradient terms
            for i = 1:element.ngrad_dofs
                dphi_i = reshape(element.cell_gradients_in_cell(edge_i, i + grad_basis_offset, iqn, :), 2, 1);
                
                % LHS
                for j = 1:element.ngrad_dofs
                    dphi_j = reshape(element.cell_gradients_in_cell(edge_i, j + grad_basis_offset, iqn, :), 2, 1);
                    MG(i,j) = MG(i,j) + wQN * (kappa_iqn * dphi_i)' * dphi_j;
                end
                
                % RHS (\GRAD vT, \GRAD w)_{PTF}
                for j = 1:element.ncell_dofs
                    dphi_j = reshape(element.cell_gradients_in_cell(edge_i, j, iqn, :), 2, 1);
                    BG(i, j) = BG(i, j) + wQN * (kappa_iqn * dphi_i)' * dphi_j;
                end
            end
            
            % Compute cell value term
            for i = 1:element.nhighorder_dofs
               phi_i = element.cell_values_in_cell(edge_i, i, iqn);
               for j = 1:element.nhighorder_dofs
                  phi_j = element.cell_values_in_cell(edge_i, j, iqn);
                  MTT(i, j) = MTT(i, j) + wQN * phi_i * phi_j;
               end
            end
        end
        
        % Compute boundary terms
        offset = element.ncell_dofs + (edge_i - 1) * element.nedge_dofs;
        
        for iqn = 1:element.nedge_qnodes
            % The quadrature point (xQN,yQN) and the weight wQN
            xQN = element.edge_quad_rules{edge_i}(iqn,1);
            yQN = element.edge_quad_rules{edge_i}(iqn,2);
            wQN = element.edge_quad_rules{edge_i}(iqn,3);

            kappa_iqn = eye(2,2);
            
            % Gradient terms
            for i = 1:element.ngrad_dofs
                dphi_i = reshape(element.cell_gradients_on_edge(edge_i, i + grad_basis_offset, iqn, :), 2, 1);
                dphi_i_nTF = dphi_i' * (kappa_iqn * nTF);
                
                % RHS (v_F, \GRAD w \dot K n_{TF})_F
                for j = 1:element.nedge_dofs
                   phi_j = element.edge_values(edge_i, j, iqn);
                   BG(i, offset + j) = BG(i, offset + j) + wQN * dphi_i_nTF * phi_j;
                end

                % RHS -(v_T, \GRAD w \dot K n_{TF})_F
                for j = 1:element.ncell_dofs
                    phi_j = element.cell_values_on_edge(edge_i, j, iqn);
                    BG(i, j) = BG(i, j) - wQN * dphi_i_nTF * phi_j;
                end
            end
            
            % Mixed (cell/edge) terms
            for i = 1:element.nedge_dofs
                phi_i = element.edge_values(edge_i, i, iqn);
                for j = 1:element.nhighorder_dofs
                   phi_j = element.cell_values_on_edge(edge_i, j, iqn);
                   MFT(edge_i,i,j) = MFT(edge_i,i,j) + wQN * phi_i * phi_j;
                end
            end
            
            % Edge/edge terms
            for i = 1:element.nedge_dofs
                phi_i = element.edge_values(edge_i, i, iqn);
                for j = 1:element.nedge_dofs
                    phi_j = element.edge_values(edge_i, j, iqn);
                    MFF(edge_i,i,j) = MFF(edge_i,i,j) + wQN * phi_i * phi_j;
                end
            end
        end    
    end
    
    % Compute the gradient reconstruction
    G_KT = MG \ BG;
    
    % Compute the stabilisation terms
    A_TF = BG' * G_KT;
    STF = zeros(element.ntotal_dofs);
    
    % High-order projection
    piRHS = MTT(1:element.ncell_dofs,(1+grad_basis_offset):element.nhighorder_dofs) * G_KT;
    piTL_rTK = MTT(1:element.ncell_dofs, 1:element.ncell_dofs) \ piRHS;
    
    for edge_i = 1:element.nedges
        hF = hho.mesh.edge_length(hho.mesh.cell_edges{cell_i}(edge_i));
        xF = hho.mesh.edge_center(hho.mesh.cell_edges{cell_i}(edge_i), :)';
        kappa_F = trace(eye(2,2));
        
        edge_terms = reshape(MFF(edge_i,:,:), element.nedge_dofs, element.nedge_dofs);
        highorder_terms = reshape(MFT(edge_i, 1:element.nedge_dofs, (1+grad_basis_offset):element.nhighorder_dofs), element.nedge_dofs, element.ngrad_dofs);
        loworder_terms = reshape(MFT(edge_i, 1:element.nedge_dofs,1:element.ncell_dofs), element.nedge_dofs, element.ncell_dofs);
        
        % High-order terms of the projection
        piF_rTK_minus_uF = edge_terms \ (highorder_terms * G_KT);
        offset = element.ncell_dofs + (edge_i-1) * element.nedge_dofs;
        piF_rTK_minus_uF(1:element.nedge_dofs, (1+offset):(offset+element.nedge_dofs)) = piF_rTK_minus_uF(1:element.nedge_dofs, (1+offset):(offset+element.nedge_dofs)) - eye(element.nedge_dofs);
        
        uT_minus_piTL_rTK = -piTL_rTK;
        uT_minus_piTL_rTK(1:element.ncell_dofs,1:element.ncell_dofs) = uT_minus_piTL_rTK(1:element.ncell_dofs,1:element.ncell_dofs) + eye(element.ncell_dofs);
        
        % Low-order terms of the projection
        piF_uT_minus_piTL_rTK = edge_terms \ (loworder_terms * uT_minus_piTL_rTK);
             
        % Stabilisation term
        BRF = piF_uT_minus_piTL_rTK + piF_rTK_minus_uF;
        STF = STF + (kappa_F / hF) * BRF' * edge_terms * BRF;
    end
    % Add stabilisation term to consistent term
    A_TF = A_TF + STF;
end

