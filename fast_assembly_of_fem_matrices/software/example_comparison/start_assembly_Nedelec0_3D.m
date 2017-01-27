
clear all
add_paths

levels=3;

load mesh3D.mat

figure(1); show_mesh(elems2nodes,nodes2coord); title('coarse mesh');

for level=0:levels
    % uniform refinement
    if (level>0)
        [nodes2coord,elems2nodes] = refinement_uniform_3D(nodes2coord,elems2nodes);
    end
    
    % get elements edgewise
    tic;
    elems2edges = get_edges(elems2nodes);
    time_edges = toc;

    %calculate affine transformations for elements
    tic;
    [B_K,~,B_K_det] = affine_transformations(nodes2coord,elems2nodes);
    time_at=toc;

    % signs in 3d for the Nedelec basis functions defined on edges
    tic
    signs = signs_edges(elems2nodes);
    time_signs = toc;
    
    % stiffness matrix assembly
    tic
    K = stiffness_matrix_Nedelec0(elems2edges,B_K,B_K_det,signs);
    time1 = toc; 
    
    % mass matrix assembly
    tic
    M = mass_matrix_Nedelec0(elems2edges,B_K,B_K_det,signs);
    time2 = toc;
    
    % print times
    rows = size(K,1);
    fprintf('time spent on level=%d, ', level);
    fprintf('edges=%f, ',time_edges);
    fprintf('affine trans=%f, ',time_at);
    fprintf('signs=%f, ',time_signs);
    fprintf('K=%f, ',time1);
    fprintf('M=%f, ',time2);
    fprintf('size of square matrices K,M=%d ',rows);
    fprintf('\n'); 
end
