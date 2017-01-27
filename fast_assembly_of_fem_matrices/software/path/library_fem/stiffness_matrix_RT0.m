function STIFF = stiffness_matrix_RT0( elems, B_K_det, signs, material )

% STIFFNESS_MATRIX_RT0
%   Vectorized calculation of stiffness matrix for the linear
%   Raviart-Thomas element in 2D/3D.
%
% SYNTAX:  STIFF = stiffness_matrix_RT0( elems, B_K_det, signs )
%          STIFF = stiffness_matrix_RT0( elems, B_K_det, signs, material )
%
% IN:   elems       elements by: edges in 2D / faces in 3D
%       B_K_det     affine map determinants
%       signs       RT basis function signs per element
%       material    elementwise constant scalar valued
%                   material coefficient
%
% OUT:  STIFF       the stiffness matrix
%

dim      = size(elems,2)-1;      % the dimension of the problem
nelems   = size(elems,1);        % number of elements
B_K_detA = abs(B_K_det);

[ip,w,nip] = intquad(1,dim);     % integration points and weighs (order 1, dimension 2/3)

% reference basis function div values on integration points,
% and number of basis functions
[~,dval,nbasis] = basis_RT0(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
STIFF = zeros(nbasis,nbasis,nelems);
for i=1:nip
    for m=1:nbasis
        for k=m:nbasis
            STIFF(m,k,:) = squeeze(STIFF(m,k,:)) + ...
                           w(i) .* B_K_detA.^(-1) .* ...
                           ( signs(:,m) .* dval(i,:,m) ) .* ...
                           ( signs(:,k) .* dval(i,:,k) );
        end
    end
end

% copy symmetric entries of the local matrices
STIFF = copy_triu(STIFF);

% apply the elemetwise constant material parameter
if ( nargin==4 )
    STIFF = astam(material,STIFF);
end

Y = reshape(repmat(elems',nbasis,1),nbasis,nbasis,nelems);    % y-indexes
X = permute(Y,[2 1 3]);                                       % x-indexes
STIFF = sparse(X(:),Y(:),STIFF(:));                           % stiffness matrix

end