function [val,cval] = exact( p )

% EXACT
%   The exact solution and its curl of the PDE: curl curl u + u = f .
%
% SYNTAX:  [val,cval] = exact( p )
%
% IN:   p         vector defining the points
%
% OUT:  val       values on points p
%       cval      curl values on points p
%

% make jumps
JUMP11 = jump('diagonal',NaN,p,'reverse');
JUMP21 = jump('diagonal',NaN,p,'reverse');

val(:,1) = JUMP11.*(sin(2.*pi.*p(:,1)) + 2.*pi.*cos(2.*pi.*p(:,1)).*(p(:,1) - p(:,2)));
val(:,2) = JUMP21.*(sin(p(:,2).*(p(:,1) - p(:,2)).^2.*(p(:,1) - 1).^2) - sin(2.*pi.*p(:,1)));

cval = 2.*JUMP11.*pi.*cos(2.*pi.*p(:,1)) - JUMP21.*(2.*pi.*cos(2.*pi.*p(:,1)) - cos(p(:,2).*(p(:,1) - p(:,2)).^2.*(p(:,1) - 1).^2).*(p(:,2).*(2.*p(:,1) - 2.*p(:,2)).*(p(:,1) - 1).^2 + p(:,2).*(2.*p(:,1) - 2).*(p(:,1) - p(:,2)).^2));