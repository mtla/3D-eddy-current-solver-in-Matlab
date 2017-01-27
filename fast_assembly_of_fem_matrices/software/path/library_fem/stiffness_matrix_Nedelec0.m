function STIFF = stiffness_matrix_Nedelec0( elems2edges, B_K, B_K_det, signs, material )

% STIFFNESS_MATRIX_NEDELEC0
%   Vectorized calculation of stiffness matrix for the 1st type
%   linear Nedelec element in 2D/3D. Note that unlike, i.e., the
%   calculation of stiffness matrix for the RT0 element, here the
%   assembly procedure differs depending on the dimension.
%
% SYNTAX:  STIFF = stiffness_matrix_Nedelec0( elems2edges, B_K, B_K_det, signs )
%          STIFF = stiffness_matrix_Nedelec0( elems2edges, B_K, B_K_det, signs, material )
%
% IN:   elems2edges    elements by edges
%       B_K            affine map matrices (used only in 3D)
%       B_K_det        affine map determinants
%       signs          Nedelec basis function signs per element
%       material       elementwise constant scalar valued
%                      material coefficient
%
% OUT:  STIFF          the stiffness matrix
%

dim      = size(B_K,1);              % the dimension of the problem
nelems   = size(elems2edges,1);      % number of elements
B_K_detA = abs(B_K_det);

[ip,w,nip] = intquad(1,dim);         % integration points and weighs (order 1, dimension 2/3)

% reference basis function curl values on integration points,
% and number of basis functions
[~,cval,nbasis] = basis_Nedelec0(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
STIFF = zeros(nbasis,nbasis,nelems);

%2d
if ( dim == 2 )
    for i=1:nip
        for m=1:nbasis
            for k=m:nbasis
                STIFF(m,k,:) = squeeze(STIFF(m,k,:)) + ...
                               w(i) .* B_K_detA.^(-1) .* ...
                               ( signs(:,m) .* cval(i,:,m) ) .* ...
                               ( signs(:,k) .* cval(i,:,k) );
            end
        end
    end

%3d
else
    for i=1:nip
        for m=1:nbasis
            for k=m:nbasis
                STIFF(m,k,:) = squeeze(STIFF(m,k,:))' + ...
                               w(i) .* B_K_detA'.^(-1) .* ...
                               sum( squeeze( astam(signs(:,m), amsv(B_K, cval(i,:,m))) ) ...
                                    .* ...
                                    squeeze( astam(signs(:,k), amsv(B_K, cval(i,:,k))) ) ...
                                  );
            end
        end
    end    
end
    
% copy symmetric entries of the local matrices
STIFF = copy_triu(STIFF);

% apply the elemetwise constant material parameter
if ( nargin==5 )
    STIFF = astam(material,STIFF);
end

Y = reshape(repmat(elems2edges',nbasis,1),nbasis,nbasis,nelems);    % y-indexes
X = permute(Y,[2 1 3]);                                             % x-indexes
STIFF = sparse(X(:),Y(:),STIFF(:));                                 % stiffness matrix

end