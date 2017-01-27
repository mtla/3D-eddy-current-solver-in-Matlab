function [val,cval,nbasis] = basis_Nedelec0(p)

% BASIS_NEDELEC0
%   The basis functions for linear Nedelec element (of first type) for
%   the space H(curl). This function evaluates these basis functions at
%   given points p in the reference geometry (triangle in 2D and
%   tetrahedron in 3D).
%
% Basis functions in 2D (edges):
%
%   N_1(x,y) = [ -y ]     N_2(x,y) = [  -y   ]     N_3(x,y) = [ 1 - y ]       
%              [  x ]                [ x - 1 ]                [   x   ]                  
%
% Basis functions in 3D (edges):
%
%   N_1(x,y,z) = [ 1 - z - y ]     N_2(x,y,z) = [     y     ]
%                [     x     ]                  [ 1 - z - x ]
%                [     x     ]                  [     y     ]
%
%   N_3(x,y,z) = [     z     ]     N_4(x,y,z) = [ - y ]
%                [     z     ]                  [   x ]
%                [ 1 - y - x ]                  [   0 ]
%
%   N_5(x,y,z) = [   0 ]           N_6(x,y,z) = [   z ]
%                [ - z ]                        [   0 ]
%                [   y ]                        [ - x ]
%
% SYNTAX:  [val,cval,nbasis] = basis_Nedelec0(p)
%
%   In the following M denotes the number of integration poins.
%
% IN:   p       M x 2/3         vector defining the points
%
% OUT:  val     M x 2/3 x 3/4   basis function values on points p
%       cval    M x 1/3 x 3/4   basis function curl values on points p
%       nbasis                  the number of basis functions
%
%                               val(j,k,i): i'th basis function
%                                           k'th component
%                                           j'th point
%

dim = size(p,2);

%2d
if ( dim == 2 )
    
    nbasis = 3;
    
    % initialize val and cval
    M    = size(p,1);
    val  = zeros( M , 2 , 3 );
    cval = zeros( M , 1 , 3 );

    % calculate basis function values
    val(:,:,1) = [ -p(:,2)    p(:,1)   ];
    val(:,:,2) = [ -p(:,2)    p(:,1)-1 ];
    val(:,:,3) = [ 1-p(:,2)   p(:,1)   ];

    % calculate basis function curl values
    cval = cval + 2;
    
%3d
elseif ( dim == 3 )
    
    nbasis = 6;
    
    % initialize val and cval
    M    = size(p,1);
    val  = zeros( M , 3 , 6 );
    cval = zeros( M , 3 , 6 );
    zero = zeros(M,1);
    two  = 2 * ones(M,1);

    % calculate basis function values
    val(:,:,1) = [ 1-p(:,3)-p(:,2)    p(:,1)             p(:,1)          ];
    val(:,:,2) = [ p(:,2)             1-p(:,3)-p(:,1)    p(:,2)          ];
    val(:,:,3) = [ p(:,3)             p(:,3)             1-p(:,2)-p(:,1) ];
    val(:,:,4) = [ -p(:,2)            p(:,1)             zero            ];
    val(:,:,5) = [ zero               -p(:,3)            p(:,2)          ];
    val(:,:,6) = [ p(:,3)             zero               -p(:,1)         ];
    
    % calculate basis function curl values
    cval(:,:,1) = [ zero    -two    two];
    cval(:,:,2) = [ two     zero    -two];
    cval(:,:,3) = [ -two    two     zero];
    cval(:,:,4) = [ zero    zero    two];
    cval(:,:,5) = [ two     zero    zero];
    cval(:,:,6) = [ zero    two     zero];
    
else
    
    error('BASIS_NEDELEC0: Input data not understood')
    
end
