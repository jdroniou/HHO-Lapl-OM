function rel_error = HHORelError( hho, u,u_exact )
%HHORelError Compute the relative error in L2 of the function
    u_interp = HHOInterpolate(hho, u_exact);
    u_diff = u - u_interp;
    abs_error = 0; % absolute error
    norm_uex = 0;
    for cell_i = 1:hho.mesh.ncells
        element = hho.elements{cell_i};
        for edge_i = 1:element.nedges
            for iqn = 1:element.ncell_qnodes
                xQN = element.cell_quad_rules{edge_i}(iqn,1);
                yQN = element.cell_quad_rules{edge_i}(iqn,2);
                wQN = element.cell_quad_rules{edge_i}(iqn,3);
                u_xQN = 0.0;
                uex_xQN = u_exact(xQN,yQN);
                for i = 1:element.ncell_dofs
                   phi_i = element.cell_values_in_cell(edge_i, i, iqn);
                   u_xQN = u_xQN + phi_i * u_diff((cell_i-1)*element.ncell_dofs+i); %u evaluated at (xQN,yQN)
                end
%                 u_diff  = u_xQN  - uex_xQN; % difference bet. u and u exact at the quadrature points
                
                abs_error = abs_error + wQN * u_xQN^2; % approximation of \int(u-u_{ex})^2
                norm_uex = norm_uex + wQN * uex_xQN^2; % approximation of \int(u_{ex})^2
            end
        end
    end
    abs_error = sqrt(abs_error);
    norm_uex = sqrt(norm_uex);
    rel_error = abs_error/norm_uex;
end

