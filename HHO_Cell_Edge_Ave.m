function [cell_val,edge_val] = HHO_Cell_Edge_Ave( hho, u )
%HHO_Cell_Edge_Ave returns 
%cell_val, the average value of the HHO function u over all cells
%edge_val, average value of u on all edges
cell_val = zeros(hho.mesh.ncells,1);
edge_val = zeros(hho.mesh.nedges,1);
ncell_dofs = hho.mesh.ncells * hho.elements{1}.ncell_dofs;

for cell_i = 1:hho.mesh.ncells
    element = hho.elements{cell_i};
    for edge_i = 1:element.nedges
        global_edge = hho.mesh.cell_edges{cell_i}(edge_i);
        if edge_val(global_edge)==0
            for s = 1:element.nedge_dofs
                for k=1:element.nedge_qnodes
                    wQN = element.edge_quad_rules{edge_i}(k,3);
                    phi_k = element.edge_values(edge_i,s,k);
                    edge_val(global_edge) = edge_val(global_edge) + wQN* u(ncell_dofs+(global_edge-1)*element.nedge_dofs+s) * phi_k;
                end
            end
            edge_val(global_edge) = edge_val(global_edge)/hho.mesh.edge_length(global_edge); % average value over an edge
        end
        
        for iqn = 1:element.ncell_qnodes
            wQN = element.cell_quad_rules{edge_i}(iqn,3);
            u_xQN = 0.0;
            for i = 1:element.ncell_dofs
                phi_i = element.cell_values_in_cell(edge_i, i, iqn);
                u_xQN = u_xQN + phi_i * u((cell_i-1)*element.ncell_dofs+i);
            end
            cell_val(cell_i) = cell_val(cell_i) + wQN * u_xQN;
        end
        
    end
    cell_val(cell_i) = cell_val(cell_i)/hho.mesh.area(cell_i); % average value over a cell
end


end

