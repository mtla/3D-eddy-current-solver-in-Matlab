fprintf('\n')
fprintf('Comparison (2D) of finite element assembly routines.\n')
fprintf('NOTE: Runtimes are reported in square brackets.\n')
fprintf('-----------------------------------------------------\n\n')

clear all
add_paths

levels = 8;

% coarse mesh of the unit square
%nodes2coord = 0.5*[0 0; 1 0; 2 0; 2 1; 1 1; 0 1; 0 2; 1 2; 2 2];
%elems2nodes = [1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7; 5 4 9; 5 9 8];
%dirichlet   = [1 2; 2 3; 3 4; 4 5; 4 9; 9 8; 8 7; 7 6; 6 1];

%switching back to Lshape as in the first paper
nodes2coord = 0.5*[0 0; 1 0; 2 0; 2 1; 1 1; 0 1; 0 2; 1 2];
elems2nodes = [1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7];
dirichlet   = [1 2];

figure(1); show_mesh(elems2nodes,nodes2coord); title('coarse mesh');

level_nodes               = zeros(levels+1,1);
level_edges               = zeros(levels+1,1);
level_time_stiffness_RT0  = zeros(levels+1,1);
level_time_mass_RT0       = zeros(levels+1,1);
level_time_stiffness_Ned0 = zeros(levels+1,1);
level_time_mass_Ned0      = zeros(levels+1,1);
level_time_stiffness_P1   = zeros(levels+1,1);
level_time_mass_P1        = zeros(levels+1,1);

for level=0:levels
    % uniform refinement
    if (level>0)
        [nodes2coord,elems2nodes,dirichlet] = refinement_uniform_2D(nodes2coord,elems2nodes,dirichlet);
    end

    % affine transformations
    tic;
    [B_K,~,B_K_det] = affine_transformations(nodes2coord,elems2nodes);
    time_at = toc;
    
    % edgewise data for Nedelec0 and RT0 element
    tic;
    elems2edges = get_edges(elems2nodes);
    signs       = signs_edges(elems2nodes);
    time_edges  = toc;
    
    % stiffness matrix assembly for RT0 element
    tic;
    K_RT0 = stiffness_matrix_RT0(elems2edges,B_K_det,signs);
    time_stiffness_RT0 = toc;
    time_stiffness_RT0 = time_stiffness_RT0 + time_at + time_edges;
    
    % mass matrix assembly for RT0 element
    tic;
    M_RT0 = mass_matrix_RT0(elems2edges,B_K,B_K_det,signs);
    time_mass_RT0 = toc;
    
    % stiffness matrix assembly for Nedelec0 element
    tic;
    K_Ned0 = stiffness_matrix_Nedelec0(elems2edges,B_K,B_K_det,signs);
    time_stiffness_Ned0 = toc;
    time_stiffness_Ned0 = time_stiffness_Ned0 + time_at + time_edges;
    
    % mass matrix assembly for Nedelec0 element
    tic;
    M_Ned0 = mass_matrix_Nedelec0(elems2edges,B_K,B_K_det,signs);
    time_mass_Ned0 = toc;
    
    % stiffness matrix assembly for P1 elements
    tic;
    K_P1 = stiffness_matrix_P1(elems2nodes,B_K,B_K_det);
    time_stiffness_P1 = toc;
    time_stiffness_P1 = time_stiffness_P1 + time_at;
    
    % mass matrix assembly for P1 elements
    tic;
    M_P1 = mass_matrix_P1(elems2nodes,B_K_det);
    time_mass_P1 = toc;
    
    level_nodes(level+1)               = max(max(elems2nodes));
    level_edges(level+1)               = max(max(elems2edges));
    level_time_stiffness_RT0(level+1)  = time_stiffness_RT0;
    level_time_mass_RT0(level+1)       = time_mass_RT0;
    level_time_stiffness_Ned0(level+1) = time_stiffness_Ned0;
    level_time_mass_Ned0(level+1)      = time_mass_Ned0;
    level_time_stiffness_P1(level+1)   = time_stiffness_P1;
    level_time_mass_P1(level+1)        = time_mass_P1;
    
    % print times
    fprintf('level=%d, '    , level);
    fprintf('elements=%d, ' , size(elems2nodes,1));
    fprintf('edges=%d, '    , max(max(elems2edges)));
    fprintf('nodes=%d \n'   , max(max(elems2nodes)));
    
    fprintf(' RT0 elements       : ');
    fprintf('time spent on K=[%f], '     , time_stiffness_RT0);
    fprintf('M=[%f], '                   , time_mass_RT0);
    fprintf(['size of matrices K,M=%d (',bytesize(K_RT0),')\n'] , size(K_RT0,1));
    
    fprintf('Ned0 elements       : ');
    fprintf('time spent on K=[%f], '     , time_stiffness_Ned0);
    fprintf('M=[%f], '                   , time_mass_Ned0);
    fprintf(['size of matrices K,M=%d (',bytesize(K_Ned0),')\n'] , size(K_Ned0,1));  
    
    fprintf('  P1 elements       : ');
    fprintf('time spent on K=[%f], '     , time_stiffness_P1);
    fprintf('M=[%f], '                   , time_mass_P1);
    fprintf(['size of matrices K,M=%d (',bytesize(K_P1),')\n'] , size(K_P1,1));
    
    fprintf('----\n');
end

% print a LaTeX-format table of the runtimes:
% level | n:o edges | RT (scale) | Nedelec (scale)

level_scale_stiffness_RT0  = zeros(levels+1,1);
level_scale_mass_RT0       = zeros(levels+1,1);
level_scale_stiffness_Ned0 = zeros(levels+1,1);
level_scale_mass_Ned0      = zeros(levels+1,1);
level_scale_stiffness_RT0(2:end,:)  = level_time_stiffness_RT0(2:end) ./ level_time_stiffness_RT0(1:end-1);
level_scale_mass_RT0(2:end,:)       = level_time_mass_RT0(2:end) ./ level_time_mass_RT0(1:end-1);
level_scale_stiffness_Ned0(2:end,:) = level_time_stiffness_Ned0(2:end) ./ level_time_stiffness_Ned0(1:end-1);
level_scale_mass_Ned0(2:end,:)      = level_time_mass_Ned0(2:end) ./ level_time_mass_Ned0(1:end-1);

for level=0:levels
    fprintf('%d ', level);
    fprintf('& ');
    fprintf('%d ', level_edges(level+1));
    fprintf('& ');
    fprintf('%2.2f (%2.2f) ', level_time_stiffness_RT0(level+1), level_scale_stiffness_RT0(level+1));
    fprintf('& ');
    fprintf('%2.2f (%2.2f) ', level_time_mass_RT0(level+1), level_scale_mass_RT0(level+1));
    fprintf('& ');
    fprintf('%2.2f (%2.2f) ', level_time_stiffness_Ned0(level+1), level_scale_stiffness_Ned0(level+1));
    fprintf('& ');
    fprintf('%2.2f (%2.2f) ', level_time_mass_Ned0(level+1), level_scale_mass_Ned0(level+1));
    fprintf('\\\\');
    fprintf('\n');
end
