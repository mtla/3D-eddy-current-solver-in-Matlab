
fprintf('\n')
fprintf('Comparison (3D) of finite element assembly routines.\n')
fprintf('NOTE: Runtimes are reported in square brackets.\n')
fprintf('-----------------------------------------------------\n\n')

clear all
add_paths

levels = 3;

load mesh3D.mat

figure(1); show_mesh(elems2nodes,nodes2coord); title('coarse mesh');

level_nodes               = zeros(levels+1,1);
level_edges               = zeros(levels+1,1);
level_faces               = zeros(levels+1,1);
level_time_stiffness_RT0  = zeros(levels+1,1);
level_time_mass_RT0       = zeros(levels+1,1);
level_time_stiffness_Ned0 = zeros(levels+1,1);
level_time_mass_Ned0      = zeros(levels+1,1);
level_time_stiffness_P1   = zeros(levels+1,1);
level_time_mass_P1        = zeros(levels+1,1);

for level=0:levels
    % uniform refinement
    if (level>0)
        [nodes2coord,elems2nodes] = refinement_uniform_3D(nodes2coord,elems2nodes);
    end
    
    % affine transformations
    tic;
    [B_K,~,B_K_det] = affine_transformations(nodes2coord,elems2nodes);
    time_at = toc;

    % facewise data for RT0 element
    tic;
    [elems2faces,faces2nodes] = get_faces(elems2nodes);
    signs_f                   = signs_faces(nodes2coord,elems2faces,faces2nodes,B_K);
    time_faces = toc;
    
    % edgewise data for Nedelec0 element
    tic;
    elems2edges = get_edges(elems2nodes);
    signs_e     = signs_edges(elems2nodes);
    time_edges  = toc;
    
    % stiffness matrix assembly for RT0 element
    tic;
    K_RT0 = stiffness_matrix_RT0(elems2faces,B_K_det,signs_f);
    time_stiffness_RT0 = toc;
    time_stiffness_RT0 = time_stiffness_RT0 + time_at + time_faces;
    
    % mass matrix assembly for RT0 element
    tic;
    M_RT0 = mass_matrix_RT0(elems2faces,B_K,B_K_det,signs_f);
    time_mass_RT0 = toc;
    
    % stiffness matrix assembly for Nedelec0 element
    tic;
    K_Ned0 = stiffness_matrix_Nedelec0(elems2edges,B_K,B_K_det,signs_e);
    time_stiffness_Ned0 = toc;
    time_stiffness_Ned0 = time_stiffness_Ned0 + time_at + time_edges;
    
    % mass matrix assembly for Nedelec0 element
    tic;
    M_Ned0 = mass_matrix_Nedelec0(elems2edges,B_K,B_K_det,signs_e);
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
    level_faces(level+1)               = max(max(elems2faces));
    level_time_stiffness_RT0(level+1)  = time_stiffness_RT0;
    level_time_mass_RT0(level+1)       = time_mass_RT0;
    level_time_stiffness_Ned0(level+1) = time_stiffness_Ned0;
    level_time_mass_Ned0(level+1)      = time_mass_Ned0;
    level_time_stiffness_P1(level+1)   = time_stiffness_P1;
    level_time_mass_P1(level+1)        = time_mass_P1;
    
    % print times
    fprintf('level=%d, '    , level);
    fprintf('elements=%d, ' , size(elems2nodes,1));
    fprintf('faces=%d, '    , max(max(elems2faces)));
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
% level | n:o faces | RT (scale) | n:o edges | Nedelec (scale)

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
    fprintf('%d ', level_faces(level+1));
    fprintf('& ');
    fprintf('%2.2f (%2.2f) ', level_time_stiffness_RT0(level+1), level_scale_stiffness_RT0(level+1));
    fprintf('& ');
    fprintf('%2.2f (%2.2f) ', level_time_mass_RT0(level+1), level_scale_mass_RT0(level+1));
    fprintf('& ');
    fprintf('%d ', level_edges(level+1));
    fprintf('& ');
    fprintf('%2.2f (%2.2f) ', level_time_stiffness_Ned0(level+1), level_scale_stiffness_Ned0(level+1));
    fprintf('& ');
    fprintf('%2.2f (%2.2f) ', level_time_mass_Ned0(level+1), level_scale_mass_Ned0(level+1));
    fprintf('\\\\');
    fprintf('\n');
end
