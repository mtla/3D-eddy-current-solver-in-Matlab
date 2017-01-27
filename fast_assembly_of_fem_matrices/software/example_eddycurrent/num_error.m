function err = num_error( elems2edges, B_K, b_K, B_K_det, signs, u, u_h )

% NUM_ERROR
%   Calculates the energy-norm of the difference between the exact solution
%   u and the numerical approximation u_h of the eddy current equation.
%
% SYNTAX:   err = num_error( elems2nodes, B_K, b_K, B_K_det, signs, u, u_h )
%
% IN:   elems2edges  elements by edges
%       B_K          affine map matrices
%       b_K          affine map vectors
%       B_K_det      affine map determinants
%       signs        signs by edges
%       u            exact solution (function handle)
%       u_h          numerical approximation coefficients
%
% OUT:  err          the difference u-u_h in energy norm (elementwise, squared)
%

dim    = size(B_K,1);
nelems = size(elems2edges,1);

B_K_detA = abs(B_K_det);

qr = 6;                        % which quadrature rule to use
[ip,w,nip] = intquad(qr,dim);  % get quadrature on the ref element

% approx. and approx. curl values
[u_h_val,u_h_cval] = solution_Nedelec0(elems2edges,B_K,B_K_det,signs,u_h,qr);

err = zeros(nelems,2);

for i=1:nip
    %calculate F_K( integration points ) in the i'th int point
    F_K_ip = squeeze(amsv(B_K, ip(i,:)))' + b_K;
    
    %get the exact solution values at the transformed integration points
    [u_val,u_cval] = u(F_K_ip);
    
    %the difference
    val_diff  = u_val  - u_h_val(:,:,i);
    cval_diff = u_cval - u_h_cval(:,:,i);
    
    err(:,1) = err(:,1) + w(i) * B_K_detA .* sum( val_diff.^2 , 2 );
    err(:,2) = err(:,2) + w(i) * B_K_detA .* sum( cval_diff.^2, 2 );
end