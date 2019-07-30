function  [grad_p]  = GradP( hho, grad_pT, cell_i, x, y)
%GradP Constructs the gradient of an HHO
%              solution to the pressure equation
%
% hho      The HHO data structure
% p        The HHO solution to the Poisson equation
% cell_i    The cell in which to compute the gradient
% (x,y)     The coordinate at which to compute the gradient
%
% The reconstruction returns the gradient of p
% as a column vector.
% 

    % Compute the gradient at the given point (x,y)
        grad = grad_pT{cell_i};
        G = [0; 0];
        for s = 1:length(grad)
            dphi_s = hho.cell_bases{cell_i}{hho.getBasisDegreeOffset(1)+s}.eval_grad(x,y);
            G = G + grad(s) * dphi_s';
        end
        grad_p = G;

end

