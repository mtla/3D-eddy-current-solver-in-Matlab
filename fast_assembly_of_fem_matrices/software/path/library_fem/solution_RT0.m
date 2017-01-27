function [sval,sdval] = solution_RT0( elems, B_K, B_K_det, signs, u_h, qr )

% SOLUTION_RT0
%   Vectorized calculation of the values of a RT0 finite element
%   approximation in the points of given integration quadrature.
%
% SYNTAX:  [sval,sdval] = solution_RT0( elems, B_K, B_K_det, signs, u_h, qr )
%
% IN:   elems     elements by: edges in 2D / faces in 3D
%       B_K       affine map matrices
%       B_K_det   affine map determinants
%       signs     RT basis function signs per element
%       u_h       RT0 finite element approximation coefficients
%       qr        quadrature rule
%
% OUT:  sval      approximation values
%       sdval     approximation divergence values
%

dim      = size(B_K,1);
nelems   = size(elems,1);

% Get integration points and basis.
% If dim==2 and qr==0, then the function was called for piecewise
% linear visualization in each element separately for a 2D function,
% so we simply take (reference) element nodes as integration points
if ( dim==2 && qr==0 )
    ip  = [0 0; 1 0; 0 1];
    nip = 3;
    [val,dval] = basis_RT0(ip);
else
    [ip,~,nip] = intquad(qr,dim);
    [val,dval] = basis_RT0(ip);
end

% initialize
sval  = zeros(nelems,dim,nip);
sdval = zeros(nelems,1  ,nip);

if ( dim==2 )
    
    for i=1:nip
        sval(:,:,i) = squeeze( ...
        	          astam(B_K_det.^(-1) .* u_h(elems(:,1)) .* signs(:,1), amsv(B_K, val(i,:,1))) + ...
                      astam(B_K_det.^(-1) .* u_h(elems(:,2)) .* signs(:,2), amsv(B_K, val(i,:,2))) + ...
                      astam(B_K_det.^(-1) .* u_h(elems(:,3)) .* signs(:,3), amsv(B_K, val(i,:,3))) )';

        sdval(:,:,i) = B_K_det.^(-1) .* ( ...
                       u_h(elems(:,1)) .* signs(:,1) .* dval(i,:,1) + ...
                       u_h(elems(:,2)) .* signs(:,2) .* dval(i,:,2) + ...
                       u_h(elems(:,3)) .* signs(:,3) .* dval(i,:,3) );
    end
    
elseif ( dim==3 )
    
    for i=1:nip
        sval(:,:,i) = squeeze( ...
        	          astam(B_K_det.^(-1) .* u_h(elems(:,1)) .* signs(:,1), amsv(B_K, val(i,:,1))) + ...
                      astam(B_K_det.^(-1) .* u_h(elems(:,2)) .* signs(:,2), amsv(B_K, val(i,:,2))) + ...
                      astam(B_K_det.^(-1) .* u_h(elems(:,3)) .* signs(:,3), amsv(B_K, val(i,:,3))) + ...
                      astam(B_K_det.^(-1) .* u_h(elems(:,4)) .* signs(:,4), amsv(B_K, val(i,:,4))) )';

        sdval(:,:,i) = B_K_det.^(-1) .* ( ...
                       u_h(elems(:,1)) .* signs(:,1) .* dval(i,:,1) + ...
                       u_h(elems(:,2)) .* signs(:,2) .* dval(i,:,2) + ...
                       u_h(elems(:,3)) .* signs(:,3) .* dval(i,:,3) + ...
                       u_h(elems(:,4)) .* signs(:,4) .* dval(i,:,4) );
    end
    
else
    
    error('SOLUTION_RT0: The input data is not understood.')
    
end