
fprintf('\n')
fprintf('Eddycurrent calculation (3D).\n')
fprintf('NOTE: Runtimes are reported in square brackets.\n')
fprintf('-----------------------------------------------\n\n')

clear all
add_paths

levels     = 0;            % how many times the mesh is refined before solving
testfolder = 'test3D_2';   % specify the test folder

% get eddycurrent problem data from the test folder
cd(testfolder)
load mesh;                 % load initial mesh
f     = @loadfun;          % right hand side
mu    = @permeability;     % permeability coefficient
kappa = @permittivity;     % permittivity coefficient

exactsolution = exist('exact','file');      % load exact solution, if it exists
if ( exactsolution ); u = @exact; end
cd ..

for i=1:levels
    [nodes2coord,elems2nodes] = refinement_uniform_3D(nodes2coord,elems2nodes);
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
[elems2faces,faces2nodes] = get_faces(elems2nodes);
bfaces2nodes              = get_boundary_faces(elems2faces,faces2nodes);
bedges2edges              = get_boundary_edges(edges2nodes,bfaces2nodes);
signs                     = signs_edges(elems2nodes);
time = toc;
fprintf('  Calculating affine transf., edge data, face data and edge signs: [%f]\n\n',time)

nelems = size(elems2nodes,1);
nnodes = size(nodes2coord,1);
nedges = max(max(elems2edges));

fprintf('  NOF elements  %d\n',nelems)
fprintf('  NOF nodes     %d\n',nnodes)
fprintf('  NOF edges     %d\n',nedges)

% solve the eddy-current problem with zero Dirichlet BC
%------------------------------------------------------
[Su,Mu,Lu] = solver_eddycurrent(elems2edges,B_K,b_K,B_K_det,signs,f,mu_const,kappa_const);
tic;
dofs       = setdiff(1:nedges,bedges2edges);
u_h        = zeros(nedges,1);
u_h(dofs)  = (Su(dofs,dofs) + Mu(dofs,dofs)) \ Lu(dofs);
time = toc;
fprintf('    Solving the linear system:  [%f]\n\n',time)

u_h_energy = 0.5*u_h'*(Su+Mu)*u_h - u_h'*Lu;
fprintf('  Energy of u_h:                  %f\n',u_h_energy)

% calculate numerical error comparing with the exact solution,
% if the exact solution was loaded
%-------------------------------------------------------------
if ( exactsolution )
    err = num_error(elems2edges,B_K,b_K,B_K_det,signs,u,u_h);
    E = sqrt(sum(sum(err)));
    fprintf('  Numerical error in energy norm:  %f\n',E)
end
fprintf('\n')
