function err = num_error( elems2nodes, B_K, b_K, B_K_det, u, u_h )

% NUM_ERROR
%   Calculates the energy-norm of the difference between the exact solution
%   u and the numerical approximation u_h of the diffusion equation.
%
% SYNTAX:   err = num_error( elems2nodes, B_K, b_K, B_K_det, u, u_h )
%
% IN:   elems2nodes  elements by nodes
%       B_K          affine map matrices
%       b_K          affine map vectors
%       B_K_det      affine map determinants
%       u            exact solution (function handle)
%       u_h          numerical approximation coefficients
%
% OUT:  err          the difference u-u_h in energy norm (elementwise, squared)
%

dim    = size(B_K,1);
nelems = size(elems2nodes,1);

B_K_detA = abs(B_K_det);

qr = 6;                        % which quadrature rule to use
[ip,w,nip] = intquad(qr,dim);  % get quadrature on the ref element

[~,u_h_gval] = solution_P1(elems2nodes,B_K,u_h,qr);  % approx. gradient values

err = zeros(nelems,1);

for i=1:nip
    %calculate F_K( integration points ) in the i'th int point
    F_K_ip = squeeze(amsv(B_K, ip(i,:)))' + b_K;
    
    %get the exact solution values at the transformed integration points
    [~,u_gval] = u(F_K_ip);
    
    %the difference
    gval_diff = u_gval - u_h_gval(:,:,i);
    
    err = err + w(i) * B_K_detA .* sum( gval_diff.^2, 2 );
end