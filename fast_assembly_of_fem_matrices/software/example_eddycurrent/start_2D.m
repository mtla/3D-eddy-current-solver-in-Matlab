
fprintf('\n')
fprintf('Eddycurrent calculation (2D).\n')
fprintf('NOTE: Runtimes are reported in square brackets.\n')
fprintf('-----------------------------------------------\n\n')

clear all
add_paths

levels     = 7;            % how many times the mesh is refined before solving
draw       = 1;            % 1=visualize; 0=no visualization
testfolder = 'test2D_3';   % specify the test folder

% get eddycurrent problem data from the test folder
cd(testfolder)
load mesh;                 % load initial mesh
load bc;                   % boundary condition
f     = @loadfun;          % right hand side
mu    = @permeability;     % permeability coefficient
kappa = @permittivity;     % permittivity coefficient

exactsolution = exist('exact','file');      % load exact solution, if it exists
if ( exactsolution ); u = @exact; end
cd ..

for i=1:levels
    [nodes2coord,elems2nodes,bedges2nodes] = refinement_uniform_2D(nodes2coord,elems2nodes,bedges2nodes);
end

% the matrix assembly routines accept only elementwise constant material
% parameters; convert the material coefficients into elementwise constant
% form if they are given as function handles
midp        = get_midpoints(elems2nodes,nodes2coord);
mu_const    = mu(midp);
kappa_const = kappa(midp);


% calculate affine transformations and sign data for Nedelec basis functions
%---------------------------------------------------------------------------
tic;
[B_K,b_K,B_K_det]         = affine_transformations(nodes2coord,elems2nodes);
[elems2edges,edges2nodes] = get_edges(elems2nodes);
bedges2edges              = get_boundary_edges(edges2nodes,bedges2nodes);
signs                     = signs_edges(elems2nodes);
time = toc;
fprintf('  Calculating affine transf., edge data, and edge signs: [%f]\n\n',time)

nelems = size(elems2nodes,1);
nnodes = size(nodes2coord,1);
nedges = max(max(elems2edges));

fprintf('  NOF elements  %d\n',nelems)
fprintf('  NOF nodes     %d\n',nnodes)
fprintf('  NOF edges     %d\n',nedges)

% solve the eddy-current problem with appropriate BC
% (full zero Dirichlet / full zero Neumann)
%---------------------------------------------------
[Su,Mu,Lu] = solver_eddycurrent(elems2edges,B_K,b_K,B_K_det,signs,f,mu_const,kappa_const);
tic;
if ( strcmp(bc,'dirichlet') )
    dofs = setdiff(1:nedges,bedges2edges);
elseif ( strcmp(bc,'neumann') )
    dofs = 1:nedges;
end
u_h        = zeros(nedges,1);
u_h(dofs)  = (Su(dofs,dofs) + Mu(dofs,dofs)) \ Lu(dofs);
time = toc;
fprintf('    Solving the linear system:  [%f]\n\n',time)

u_h_energy = 0.5*u_h'*(Su+Mu)*u_h - u_h'*Lu;
fprintf('  Energy of u_h:                  %f\n',u_h_energy)

% calculate numerical error comparing with the exact solution
%------------------------------------------------------------
if ( exactsolution )
    err = num_error(elems2edges,B_K,b_K,B_K_det,signs,u,u_h);
    E = sqrt(sum(sum(err)));
    fprintf('  Numerical error in energy norm:  %f\n',E)
end
fprintf('\n')

% visualize functions
%--------------------
if ( draw == 1 )
    % approximation u_h
    % note that the approximation is elementwise linear but the visualization
    % is done by taking a constant value (the midpoint) on each element
    %------------------------------------------------------------------------
    
    [u_h_val,u_h_cval] = solution_Nedelec0(elems2edges,B_K,B_K_det,signs,u_h,1);
    u_h_val_x = u_h_val(:,1,:);
    u_h_val_y = u_h_val(:,2,:);
    
    figure
    show_mesh(elems2nodes,nodes2coord)
    
    figure
    show_constant_scalar(u_h_val_x,nodes2coord,elems2nodes)
    axis on
    title('edge approximation of the eddy current problem (x-component)')
    view([-73,40]);

    figure
    show_constant_scalar(u_h_val_y,nodes2coord,elems2nodes)
    axis on
    title('edge approximation of the eddy current problem (y-component)')
    view([-63,42]);
    
    figure
    show_constant_scalar(u_h_cval,nodes2coord,elems2nodes)
    axis on
    title('curl of edge approximation of the eddy current problem')
    view([-49,50]);
end
