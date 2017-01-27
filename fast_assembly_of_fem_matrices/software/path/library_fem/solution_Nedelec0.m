function [sval,scval] = solution_Nedelec0( elems2edges, B_K, B_K_det, signs, u_h, qr )

% SOLUTION_Nedelec0
%   Vectorized calculation of the values of a Nedelec0 finite element
%   approximation in the points of given integration quadrature.
%
% SYNTAX:  [sval,scval] = solution_Nedelec0( elems2edges, B_K, B_K_det, signs, u_h, qr )
%
% IN:   elems     elements by edges (both in 2D/3D)
%       B_K       affine map matrices
%       B_K_det   affine map determinants
%       signs     Nedelec basis function signs per element
%       u_h       Nedelec0 finite element approximation coefficients
%       qr        quadrature rule
%
% OUT:  sval      approximation values
%       scval     approximation curl values
%

dim      = size(B_K,1);
nelems   = size(elems2edges,1);
B_K_invT = amt(aminv(B_K));

% Get integration points and basis.
% If dim==2 and qr==0, then the function was called for piecewise
% linear visualization in each element separately for a 2D function,
% so we simply take (reference) element nodes as integration points
if ( dim==2 && qr==0 )
    ip  = [0 0; 1 0; 0 1];
    nip = 3;
    [val,cval] = basis_Nedelec0(ip);
else
    [ip,~,nip] = intquad(qr,dim);
    [val,cval] = basis_Nedelec0(ip);
end

if ( dim==2 )
    % initialize
    sval  = zeros(nelems,dim,nip);
    scval = zeros(nelems,1  ,nip);
    
    for i=1:nip
        sval(:,:,i) = squeeze( ...
        	          astam(u_h(elems2edges(:,1)) .* signs(:,1), amsv(B_K_invT, val(i,:,1))) + ...
                      astam(u_h(elems2edges(:,2)) .* signs(:,2), amsv(B_K_invT, val(i,:,2))) + ...
                      astam(u_h(elems2edges(:,3)) .* signs(:,3), amsv(B_K_invT, val(i,:,3))) )';

        scval(:,:,i) = B_K_det.^(-1) .* ( ...
                       u_h(elems2edges(:,1)) .* signs(:,1) .* cval(i,:,1) + ...
                       u_h(elems2edges(:,2)) .* signs(:,2) .* cval(i,:,2) + ...
                       u_h(elems2edges(:,3)) .* signs(:,3) .* cval(i,:,3) );
    end
    
elseif ( dim==3 )
    % initialize
    sval  = zeros(nelems,dim,nip);
    scval = zeros(nelems,dim,nip);
    
    for i=1:nip
        sval(:,:,i) = squeeze( ...
        	          astam(u_h(elems2edges(:,1)) .* signs(:,1), amsv(B_K_invT, val(i,:,1))) + ...
                      astam(u_h(elems2edges(:,2)) .* signs(:,2), amsv(B_K_invT, val(i,:,2))) + ...
                      astam(u_h(elems2edges(:,3)) .* signs(:,3), amsv(B_K_invT, val(i,:,3))) + ...
                      astam(u_h(elems2edges(:,4)) .* signs(:,4), amsv(B_K_invT, val(i,:,4))) + ...
                      astam(u_h(elems2edges(:,5)) .* signs(:,5), amsv(B_K_invT, val(i,:,5))) + ...
                      astam(u_h(elems2edges(:,6)) .* signs(:,6), amsv(B_K_invT, val(i,:,6))) )';

        scval(:,:,i) = squeeze( ...
                       astam(B_K_det.^(-1) .* u_h(elems2edges(:,1)) .* signs(:,1), amsv(B_K, cval(i,:,1))) + ...
                       astam(B_K_det.^(-1) .* u_h(elems2edges(:,2)) .* signs(:,2), amsv(B_K, cval(i,:,2))) + ...
                       astam(B_K_det.^(-1) .* u_h(elems2edges(:,3)) .* signs(:,3), amsv(B_K, cval(i,:,3))) + ...
                       astam(B_K_det.^(-1) .* u_h(elems2edges(:,4)) .* signs(:,4), amsv(B_K, cval(i,:,4))) + ...
                       astam(B_K_det.^(-1) .* u_h(elems2edges(:,5)) .* signs(:,5), amsv(B_K, cval(i,:,5))) + ...
                       astam(B_K_det.^(-1) .* u_h(elems2edges(:,6)) .* signs(:,6), amsv(B_K, cval(i,:,6))) )';
    end
    
else
    
    error('SOLUTION_NEDELEC0: The input data is not understood.')
    
end