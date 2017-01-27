clear all
add_paths

levels=6;

% coarse mesh of the unit square
nodes2coord = 0.5*[0 0; 1 0; 2 0; 2 1; 1 1; 0 1; 0 2; 1 2; 2 2];
elems2nodes = [1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7; 5 4 9; 5 9 8];
dirichlet   = [1 2; 2 3; 3 4; 4 9; 9 8; 8 7; 7 6; 6 1];

figure(1); show_mesh(elems2nodes,nodes2coord); title('coarse mesh');

for level=0:levels
    % uniform refinement
    if (level>0)
        [nodes2coord,elems2nodes,dirichlet] = refinement_uniform_2D(nodes2coord,elems2nodes,dirichlet);
    end
    
    %calculate affine transformations for elements
    tic;
    [B_K,~,B_K_det] = affine_transformations(nodes2coord,elems2nodes);
    time_at = toc;
    
    % stiffness matrix assembly
    tic
    K = stiffness_matrix_P1(elems2nodes,B_K,B_K_det);
    time1 = toc;
    
    % mass matrix assembly
    tic;
    M = mass_matrix_P1(elems2nodes,B_K_det);
    time2 = toc;
    
    % print times 
    rows = size(K,1);
    fprintf('time spent on level=%d, ', level);
    fprintf('affine trans=%f, ',time_at);
    fprintf('K=%f, ',time1);
    fprintf('M=%f, ',time2);
    fprintf('size of square matrices K,M=%d ',rows);
    fprintf('\n');
end
