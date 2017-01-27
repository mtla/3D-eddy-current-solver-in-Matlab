function MASS = mass_matrix_Nedelec0( elems2edges, B_K, B_K_det, signs, material )

% MASS_MATRIX_NEDELEC0
%   Vectorized calculation of mass matrix for the 1st type
%   linear Nedelec element in 2D/3D.
%
% SYNTAX:  MASS = mass_matrix_Nedelec0( elems2edges, B_K, B_K_det, signs )
%          MASS = mass_matrix_Nedelec0( elems2edges, B_K, B_K_det, signs, material )
%
% IN:   elems2edges    elements by edges (in both 2D/3D)
%       B_K            affine map matrices
%       B_K_det        affine map determinants
%       signs          Nedelec basis function signs per element
%       material       elementwise constant scalar valued
%                      material coefficient
%
% OUT:  MASS           the mass matrix
%

dim      = size(B_K,1);              % the dimension of the problem
nelems   = size(elems2edges,1);      % number of elements
B_K_detA = abs(B_K_det);
B_K_invT = amt(aminv(B_K));

[ip,w,nip] = intquad(2,dim);         % integration points and weighs (order 2, dimension 2/3)

% reference basis function curl values on integration points,
% and number of basis functions
[val,~,nbasis] = basis_Nedelec0(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
MASS = zeros(nbasis,nbasis,nelems);
for i=1:nip
    for m=1:nbasis
        for k=m:nbasis
            MASS(m,k,:) = squeeze(MASS(m,k,:))' + ...
                          w(i) .* B_K_detA' .* ...
                          sum( squeeze( astam(signs(:,m), amsv(B_K_invT, val(i,:,m))) ) ...
                               .* ...
                               squeeze( astam(signs(:,k), amsv(B_K_invT, val(i,:,k))) ) ...
                             );
        end
    end
end

% copy symmetric entries of the local matrices
MASS = copy_triu(MASS);

% apply the elemetwise constant material parameter
if ( nargin==5 )
    MASS = astam(material,MASS);
end

Y = reshape(repmat(elems2edges',nbasis,1),nbasis,nbasis,nelems);    % y-indexes
X = permute(Y,[2 1 3]);                                             % x-indexes
MASS = sparse(X(:),Y(:),MASS(:));                                   % mass matrix

end