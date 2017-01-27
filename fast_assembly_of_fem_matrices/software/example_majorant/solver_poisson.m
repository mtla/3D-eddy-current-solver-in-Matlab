function [STIFF,LOAD] = solver_poisson( elems2nodes, B_K, b_K, B_K_det, f )

% SOLVER_POISSON
%   Linear H^1-FEM solver in 2D / 3D.
%
% SYNTAX:  solver(elems2nodes, B_K, b_K, B_K_det, f)
%
% IN:   elems2nodes  elements by nodes (in both 2D/3D)
%       B_K          affine map matrices
%       b_K          affine map vectors
%       B_K_det      affine map determinants
%       f            load function
%
% OUT:  STIFF        stiffness matrix
%       LOAD         load vector
%


% calculate the local STIFFNESS matrices
%---------------------------------------
tic;
STIFF = stiffness_matrix_P1( elems2nodes, B_K, B_K_det );
time_STIFF = toc;


%calculate the local LOAD vectors
%--------------------------------
tic;

dim    = size(B_K,1);
nelems = size(elems2nodes,1);

B_K_detA = abs(B_K_det);

% get the integration points and weighs on the ref element
[ipL,wL,nipL] = intquad(3,dim);        % L as in load

% get the ref basis function and grad values on the int points
[valL,~,nbasis] = basis_P1(ipL);

LOAD = zeros(nelems,nbasis);

for i=1:nipL
    %calculate F_K( integration points ) in the i'th int point
    F_K_ip = squeeze(amsv(B_K, ipL(i,:)))' + b_K;
    
    %get the load function values at these transformed integration points
    fval = f(F_K_ip);
    
    for k=1:nbasis
        LOAD(:,k) = LOAD(:,k) + ...
            	    wL(i) .* B_K_detA .* ...
                    ( fval .* valL(i,:,k) );
    end
end
LOAD = sparse(elems2nodes,1,LOAD);      % assemble data into full vectors

time_LOAD = toc;


% print runtimes
%---------------
fprintf('\n')
fprintf('  SOLVER_POISSON\n')
fprintf('    STIFFNESS matrix:       [%f]\n',time_STIFF)
fprintf('    LOAD vector:            [%f]\n',time_LOAD )
fprintf('\n')

end