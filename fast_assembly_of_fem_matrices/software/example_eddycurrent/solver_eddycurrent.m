function [STIFF,MASS,LOAD] = solver_eddycurrent(elems2edges, B_K, b_K, B_K_det, signs, f, mu, kappa)

% SOLVER_EDDYCURRENT
%   Linear Nedelec-FEM solver in 2D / 3D.
%
% SYNTAX:   [STIFF,MASS,LOAD] = solver(elems2edges, B_K, b_K, B_K_det, signs, f, mu, kappa)
%
% IN:   elems2edges  elements by edges (in both 2D/3D)
%       B_K          affine map matrices
%       b_K          affine map vectors
%       B_K_det      affine map determinants
%       signs        edge signs
%       f            load function
%       mu           elementwise constant scalar valued permeability coefficient
%       kappa        elementwise constant scalar valued permittivity coefficient
%
% OUT:  STIFF        stiffness matrix
%       MASS         mass matrix
%       LOAD         load vector
%

dim = size(B_K,1);


% calculate the local STIFFNESS and MASS matrices
%------------------------------------------------
tic;
STIFF = stiffness_matrix_Nedelec0(elems2edges,B_K,B_K_det,signs,mu.^(-1));
time_STIFF = toc; 

tic;
MASS = mass_matrix_Nedelec0(elems2edges,B_K,B_K_det,signs,kappa);
time_MASS = toc;


%calculate the local LOAD vectors
%--------------------------------
tic;

nelems   = size(elems2edges,1);
B_K_detA = abs(B_K_det);
B_K_invT = amt(aminv(B_K));

% get the integration points and weighs on the ref element
[ipL,wL,nipL] = intquad(6,dim);        % L as in load

% get the ref basis function and grad values on the int points
[valL,~,nbasis] = basis_Nedelec0(ipL);

LOAD = zeros(nelems,nbasis);

for i=1:nipL
    %calculate F_K( integration points ) in the i'th int point
    F_K_ip = squeeze(amsv(B_K, ipL(i,:)))' + b_K;
    
    %get the load function values at these transformed integration points
    fval = f(F_K_ip);
    
    for k=1:nbasis
        LOAD(:,k) = LOAD(:,k) + ...
                    wL(i) .* B_K_detA .* ...
                    sum( fval ...
                         .* ...
                         squeeze( astam(signs(:,k), amsv(B_K_invT, valL(i,:,k))) )' ...
                       , 2);
    end
end
LOAD = sparse(elems2edges,1,LOAD);        % assemble data into full vectors

time_LOAD = toc;


% print runtimes
%---------------
fprintf('\n')
fprintf('  SOLVER_EDDYCURRENT\n')
fprintf('    STIFFNESS matrix:       [%f]\n',time_STIFF)
fprintf('    MASS matrix:            [%f]\n',time_MASS)
fprintf('    LOAD vector:            [%f]\n',time_LOAD)
fprintf('\n')

end