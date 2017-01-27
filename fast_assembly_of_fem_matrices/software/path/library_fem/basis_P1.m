function [val,gval,nbasis] = basis_P1(p)

% BASIS_P1
%   The basis functions for linear Courant/Lagrange/Nodal element for
%   the space H^1. This function evaluates these basis functions at
%   given points p in the reference geometry (triangle in 2D and
%   tetrahedron in 3D).
%
% Basis functions in 2D (nodes):
%
%   N_1(x,y) = 1 - x - y        N_2(x,y) = x        N_3(x,y) = y
%
% Basis functions in 3D (nodes):
%
%   N_1(x,y,z) = 1 - x - y - z          N_2(x,y,z) = x
%   N_3(x,y,z) = y                      N_4(x,y,z) = z
%
% SYNTAX:  [val,gval,nbasis] = basis_P1(p)
%
%   In the following M denotes the number of integration points in p.
%
% IN:   p       M x 2/3         vector defining the points
%
% OUT:  val     M x  1  x 3/4   basis function values on points p
%       gval    M x 2/3 x 3/4   basis function grad values on points p
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

    % initialize val and gval tensor
    M = size(p,1);
    val  = zeros( M , 1 , nbasis );
    gval = zeros( M , 2 , nbasis );
    unit = ones(M,1);          %help var
    zero = zeros(M,1);

    % calculate basis function values
    val(:,:,1) = -p(:,1)-p(:,2)+1 ;
    val(:,:,2) = p(:,1) ;
    val(:,:,3) = p(:,2) ;

    % calculate gradient(basis function) values
    gval(:,:,1) = [ -unit  -unit ];
    gval(:,:,2) = [  unit   zero ];
    gval(:,:,3) = [  zero   unit ];
    
%3d
elseif ( dim == 3 )
    
    nbasis = 4;

    % initialize val and gval tensor
    M = size(p,1);
    val  = zeros( M , 1 , nbasis );
    gval = zeros( M , 3 , nbasis );
    unit = ones(M,1);
    zero = zeros(M,1);

    % calculate basis function values
    val(:,:,1) = 1-p(:,1)-p(:,2)-p(:,3) ;
    val(:,:,2) = p(:,1) ;
    val(:,:,3) = p(:,2) ;
    val(:,:,4) = p(:,3) ;

    % calculate gradient(basis function) values
    gval(:,:,1) = [ -unit  -unit  -unit ];
    gval(:,:,2) = [  unit   zero   zero ];
    gval(:,:,3) = [  zero   unit   zero ];
    gval(:,:,4) = [  zero   zero   unit ];
    
else
    
    error('BASIS_P1: Input data not understood')
    
end
