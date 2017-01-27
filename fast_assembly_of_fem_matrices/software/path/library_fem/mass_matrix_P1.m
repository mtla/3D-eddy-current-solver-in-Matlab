function MASS = mass_matrix_P1( elems2nodes, B_K_det, material )

% MASS_MATRIX_P1
%   Vectorized calculation of mass matrix for the linear Courant (P1)
%   element in 2D/3D.
%
% SYNTAX:  MASS = mass_matrix_P1( elems2nodes, B_K_det )
%          MASS = mass_matrix_P1( elems2nodes, B_K_det, material )
%
% IN:   elems2nodes    elements by nodes (in both 2D/3D)
%       B_K_det        affine map determinants
%       material       elementwise constant scalar valued
%                      material coefficient
%
% OUT:  MASS           the mass matrix
%

dim      = size(elems2nodes,2)-1;      % the dimension of the problem
nelems   = size(elems2nodes,1);        % number of elements
B_K_detA = abs(B_K_det);

[ip,w,nip] = intquad(2,dim);     % integration points and weighs (order 2, dimension 2/3)

% reference basis function values on integration points,
% and number of basis functions
[val,~,nbasis] = basis_P1(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
MASS = zeros(nbasis,nbasis,nelems);
for i=1:nip
    for m=1:nbasis
        for k=m:nbasis   
            MASS(m,k,:) = squeeze(MASS(m,k,:)) + ...
                          w(i) .* B_K_detA .* ...
                          ( val(i,:,m) .* val(i,:,k) );
        end
    end
end

% copy symmetric entries of the local matrices
MASS = copy_triu(MASS);

% apply the elemetwise constant material parameter
if ( nargin==3 )
    MASS = astam(material,MASS);
end

Y = reshape(repmat(elems2nodes',nbasis,1),nbasis,nbasis,nelems);  % y-indexes
X = permute(Y,[2 1 3]);                                           % x-indexes
MASS = sparse(X(:),Y(:),MASS(:));                                 % mass matrix

end