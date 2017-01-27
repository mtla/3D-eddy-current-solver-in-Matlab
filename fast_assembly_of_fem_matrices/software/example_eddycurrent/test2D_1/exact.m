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

val(:,1) = p(:,2).*(1 - p(:,2));
val(:,2) = p(:,1).*(1 - p(:,1));

cval = 2.*p(:,2) - 2.*p(:,1);