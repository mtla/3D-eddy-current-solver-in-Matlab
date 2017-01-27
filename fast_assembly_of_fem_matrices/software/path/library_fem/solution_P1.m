function [sval,sgval] = solution_P1( elems2nodes, B_K, u_h, qr)

% SOLUTION_P1
%   Vectorized calculation of the values of a P1 finite element
%   approximation in the points of given integration quadrature.
%
% SYNTAX:  [sval,sgval] = solution_P1( elems2nodes, B_K, u_h, qr)
%
% IN:   elems2nodes    elements by nodes (in both 2D/3D)
%       B_K            affine map matrices
%       u_h            P1 finite element approximation coefficients
%       qr             quadrature rule
%
% OUT:  sval           approximation values
%       sgval          approximation gradient values
%

dim      = size(B_K,1);
nelems   = size(elems2nodes,1);
B_K_invT = amt(aminv(B_K));

% Get integration points and basis.
% If dim==2 and qr==0, then the function was called for piecewise
% linear visualization in each element separately for a 2D function,
% so we simply take (reference) element nodes as integration points
if ( dim==2 && qr==0 )
    ip  = [0 0; 1 0; 0 1];
    nip = 3;
    [val,gval] = basis_P1(ip);
else
    [ip,~,nip] = intquad(qr,dim);
    [val,gval] = basis_P1(ip);
end

% initialize
sval  = zeros(nelems,1  ,nip);
sgval = zeros(nelems,dim,nip);

if ( dim==2 )
    
    for i=1:nip
        sval(:,:,i) = u_h(elems2nodes(:,1)) .* val(i,:,1) + ...
                      u_h(elems2nodes(:,2)) .* val(i,:,2) + ...
                      u_h(elems2nodes(:,3)) .* val(i,:,3) ;

        sgval(:,:,i) = squeeze( ...
                       astam(u_h(elems2nodes(:,1)), amsv(B_K_invT, gval(i,:,1))) + ...
                       astam(u_h(elems2nodes(:,2)), amsv(B_K_invT, gval(i,:,2))) + ...
                       astam(u_h(elems2nodes(:,3)), amsv(B_K_invT, gval(i,:,3))) )';
    end
    
elseif ( dim==3 )
    
    for i=1:nip
        sval(:,:,i) = u_h(elems2nodes(:,1)) .* val(i,:,1) + ...
                      u_h(elems2nodes(:,2)) .* val(i,:,2) + ...
                      u_h(elems2nodes(:,3)) .* val(i,:,3) + ...
                      u_h(elems2nodes(:,4)) .* val(i,:,4);

        sgval(:,:,i) = squeeze( ...
                       astam(u_h(elems2nodes(:,1)), amsv(B_K_invT, gval(i,:,1))) + ...
                       astam(u_h(elems2nodes(:,2)), amsv(B_K_invT, gval(i,:,2))) + ...
                       astam(u_h(elems2nodes(:,3)), amsv(B_K_invT, gval(i,:,3))) + ...
                       astam(u_h(elems2nodes(:,4)), amsv(B_K_invT, gval(i,:,4))) )';
    end
    
else
    
    error('SOLUTION_P1: The input data is not understood.')
    
end