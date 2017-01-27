function [val,gval] = exact( p )

% EXACT
%   The exact solution and its grad of the PDE: - div grad u = f .
%
% SYNTAX:  [val,gval] = exact( p )
%
% IN:   p         vector defining the points
%
% OUT:  val       values on points p
%       cval      grad values on points p
%

val = p(:,1).*p(:,2).*(p(:,1) - 1).*(p(:,2) - 1);

gval = [ ...
p(:,2).*(2.*p(:,1) - 1).*(p(:,2) - 1) ...
p(:,1).*(2.*p(:,2) - 1).*(p(:,1) - 1) ...
];