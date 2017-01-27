function [STIFF,MASS,LOAD1,LOAD2] = solver_majorant(elems2nodes, elems, B_K, b_K, B_K_det, signs, f, u_h)

% SOLVER_MAJORANT
%   Linear RT0-FEM solver in 2D / 3D.
%
% SYNTAX:  solver(elems2nodes, elems, B_K, b_K, B_K_det, signs, f, u_h)
%
% IN:   elems2nodes  elements by nodes (in both 2D/3D)
%       elems        elements by: edges in 2D / faces in 3D
%       B_K          affine map matrices
%       b_K          affine map vectors
%       B_K_det      affine map determinants
%       signs        RT0 basis function signs per element
%       f            load function
%       u_h          poisson problem approximation coefficients
%
% OUT:  STIFF        stiffness matrix
%       MASS         mass matrix
%       LOAD1        load vector part 1
%       LOAD2        load vector part 2
%


% calculate the local STIFFNESS and MASS matrices
%------------------------------------------------
tic;
STIFF = stiffness_matrix_RT0( elems, B_K_det, signs );
time_STIFF = toc;

tic;
MASS = mass_matrix_RT0( elems, B_K, B_K_det, signs );
time_MASS = toc;


% calculate the local LOAD vectors
%---------------------------------
tic;

dim    = size(B_K,1);
nelems = size(elems,1);

B_K_detA = abs(B_K_det);

% get the integration points and weighs on the ref element
[ipD,wD,nipD] = intquad(1,dim);             % D as in differentials
[ipL,wL,nipL] = intquad(6,dim);             % L as in load

% get the ref basis function and their div values on the integration points
[valD,~    ,nbasis] = basis_RT0(ipD);
[~   ,dvalL,~     ] = basis_RT0(ipL);

% initialize matrices and vectors
LOAD1 = zeros(nelems,nbasis);
LOAD2 = LOAD1;

for i=1:nipL
    %calculate F_K( integration points ) in the i'th int point
    F_K_ip = squeeze(amsv(B_K, ipL(i,:)))' + b_K;
    
    %get the load function values at these transformed integration points
    fval = f(F_K_ip);
    
    for k=1:nbasis
        LOAD1(:,k) = LOAD1(:,k) + ...
                     wL(i) .* B_K_detA .* ...
                     fval .* signs(:,k) .* B_K_det.^(-1) .* dvalL(i,:,k);
    end    
end
LOAD1 = sparse(elems,1,LOAD1);  % assemble data into full vectors

[~,u_h_gval] = solution_P1(elems2nodes,B_K,u_h,1);  % solution gradient values
for i=1:nipD
    for k=1:nbasis
        LOAD2(:,k) = LOAD2(:,k) + ...
                     wD(i) .* B_K_detA .* ...
                     sum( u_h_gval(:,:,i) .* ...
                          squeeze(astam(signs(:,k).*B_K_det.^(-1), amsv(B_K, valD(i,:,k))))' ...
                        , 2);
    end
end
LOAD2 = sparse(elems,1,LOAD2);  % assemble data into full vectors

time_LOAD = toc;


% print runtimes
%---------------
fprintf('\n')
fprintf('  SOLVER_MAJORANT\n')
fprintf('    STIFFNESS matrix:       [%f]\n',time_STIFF)
fprintf('    MASS matrix:            [%f]\n',time_MASS )
fprintf('    LOAD vector:            [%f]\n',time_LOAD )
fprintf('\n')

end