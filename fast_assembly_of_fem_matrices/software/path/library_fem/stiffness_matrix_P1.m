function STIFF = stiffness_matrix_P1( elems2nodes, B_K, B_K_det, material )

% STIFFNESS_MATRIX_P1
%   Vectorized calculation of stiffness matrix for the linear Courant (P1)
%   element in 2D/3D.
%
% SYNTAX:  STIFF = stiffness_matrix_P1( elems2nodes, B_K, B_K_det )
%          STIFF = stiffness_matrix_P1( elems2nodes, B_K, B_K_det, material )
%
% IN:   elems2nodes    elements by nodes (in both 2D/3D)
%       B_K            affine map matrices
%       B_K_det        affine map determinants
%       material       elementwise constant scalar valued
%                      material coefficient
%
% OUT:  STIFF          the stiffness matrix
%

dim      = size(B_K,1);          % the dimension of the problem
nelems   = size(elems2nodes,1);  % number of elements
B_K_detA = abs(B_K_det);
B_K_invT = amt(aminv(B_K));

[ip,w,nip] = intquad(1,dim);     % integration points and weighs (order 1, dimension 2/3)

% reference basis function grad values on integration points,
% and number of basis functions
[~,gval,nbasis] = basis_P1(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
STIFF = zeros(nbasis,nbasis,nelems);
for i=1:nip
    for m=1:nbasis
        for k=m:nbasis 
            STIFF(m,k,:) = squeeze(STIFF(m,k,:))' + ...
                           w(i) .* B_K_detA' .* ...
                           sum( squeeze(amsv(B_K_invT, gval(i,:,m))) ...
                                .* ...
                                squeeze(amsv(B_K_invT, gval(i,:,k))) ...
                              );
        end
    end
end

% copy symmetric entries of the local matrices
STIFF = copy_triu(STIFF);

% apply the elemetwise constant material parameter
if ( nargin==4 )
    STIFF = astam(material,STIFF);
end

Y = reshape(repmat(elems2nodes',nbasis,1),nbasis,nbasis,nelems);  % y-indexes
X = permute(Y,[2 1 3]);                                           % x-indexes
STIFF = sparse(X(:),Y(:),STIFF(:));                               % stiffness matrix

end