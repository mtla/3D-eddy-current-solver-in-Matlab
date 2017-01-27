
clear all
add_paths

levels=6;

% coarse mesh of the unit square
coordinates = 0.5*[0 0; 1 0; 2 0; 2 1; 1 1; 0 1; 0 2; 1 2; 2 2];
elems2nodes = [1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7; 5 4 9; 5 9 8];
dirichlet   = [1 2; 2 3; 3 4; 4 9; 9 8; 8 7; 7 6; 6 1];

figure(1); show_mesh(elems2nodes,coordinates); title('coarse mesh');

for level=0:levels
    % uniform refinement
    if (level>0)
        [coordinates,elems2nodes,dirichlet] = refinement_uniform_2D(coordinates,elems2nodes,dirichlet);
    end

    % get elements edgewise
    tic;
    elems2edges = get_edges(elems2nodes);
    time_edges = toc;

    %calculate affine transformations for elements
    tic;
    [B_K,~,B_K_det] = affine_transformations(coordinates,elems2nodes);
    time_at = toc;

    % signs in 2d for the Nedelec basis functions defined on edges
    tic;
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
