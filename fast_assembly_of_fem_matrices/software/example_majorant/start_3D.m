
fprintf('\n')
fprintf('Majorant calculation (3D) for a given approximation of the Poisson problem.\n')
fprintf('NOTE: Runtimes are reported in square brackets.\n')
fprintf('---------------------------------------------------------------------------\n\n')

clear all
add_paths

level      = 0;            % how many times the mesh is refined before solving
testfolder = 'test3D';     % specify the test folder

% get Poisson problem data from the test folder
cd(testfolder)
load mesh;                 % load initial mesh
f = @loadfun;              % right hand side
u = @exact;                % exact solution
load const;                % load Cf, the Friedrichs constant
cd ..

for i=1:level
    [nodes2coord,elems2nodes] = refinement_uniform_3D(nodes2coord,elems2nodes);
end

% calculate affine transformations and sign data for RT basis functions
%----------------------------------------------------------------------
tic;
[B_K,b_K,B_K_det]         = affine_transformations(nodes2coord,elems2nodes);
[elems2faces,faces2nodes] = get_faces(elems2nodes);
bfaces2nodes              = get_boundary_faces(elems2faces,faces2nodes);
signs                     = signs_faces(nodes2coord,elems2faces,faces2nodes,B_K);
time = toc;
fprintf('  Calculating affine transf., face data, and face signs: [%f]\n\n',time)

nelems = size(elems2nodes,1);
nnodes = size(nodes2coord,1);
nfaces = size(faces2nodes,1);

fprintf('  NOF elements  %d\n',nelems)
fprintf('  NOF nodes     %d\n',nnodes)
fprintf('  NOF faces     %d\n',nfaces)

% solve the Poisson problem with zero Dirichlet BC
%-------------------------------------------------
[Su,Lu] = solver_poisson(elems2nodes,B_K,b_K,B_K_det,f);
tic;
bnodes2nodes    = unique(bfaces2nodes);
dofs            = setdiff(1:nnodes,bnodes2nodes);
u_h             = zeros(nnodes,1);
u_h(dofs)       = Su(dofs,dofs) \ Lu(dofs);
time = toc;
fprintf('    Solving the linear system:  [%f]\n\n',time)

% calculate numerical error comparing with the exact solution
err = num_error( elems2nodes, B_K, b_K, B_K_det, u, u_h );
E = sqrt(sum(err));
fprintf('  Numerical error in energy norm: %f\n\n',E)

% calculate the majorant
%-----------------------
[Sy,My,L1y,L2y] = solver_majorant(elems2nodes, elems2faces, B_K, b_K, B_K_det, signs, f, u_h);

tic;
% these are needed for calculating the majorant value with matrix operations
uhSuh = u_h'*Su*u_h;
f_L2norm = norm_L2(elems2nodes,B_K,b_K,B_K_det,f);
beta = 1;
for i=1:10
    % solve the matrix system
    y_h = (  (1+beta) * Cf^2 * Sy  + (1+1/beta) * My  ) \ ...
          ( -(1+beta) * Cf^2 * L1y + (1+1/beta) * L2y );
    y_h = full(y_h);
          
	% calculate the two terms of the majorant with matrix operations
    maj_residual = f_L2norm + y_h'*Sy*y_h + 2*y_h'*L1y;
    maj_dual     = uhSuh + y_h'*My*y_h - 2*y_h'*L2y;
    
    % calculate the two terms of the majorant by integration
    %[maj_residual,maj_dual] = majorant(elems2nodes, elems2faces, B_K, b_K, B_K_det, signs, f, u_h, y_h);
    
    maj = sqrt( (1+beta)*Cf^2*sum(maj_residual) + (1+1/beta)*sum(maj_dual) );
    
    % efficiency index of the majorant
    Ieff = maj/E;
    
    % print
    fprintf('    Iteration %d ',i);
    fprintf(': beta=%f ' ,beta);
    fprintf('; maj=%f '  ,maj );
    fprintf('; Ieff=%f\n',Ieff);
    
    % stopping criterion
    if ( i>1 && (maj_prev-maj)/maj_prev < 10^(-4) )
        break;
    end
    
    beta = sqrt(sum(maj_dual))/(Cf*sqrt(sum(maj_residual)));
    maj_prev = maj;
end
time = toc;
fprintf('    Time taken for iterations:  [%f]\n\n',time)
