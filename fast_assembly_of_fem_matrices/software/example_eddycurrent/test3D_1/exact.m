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

val  = zeros(size(p));
cval = zeros(size(p));

val(:,3) = sin(pi.*p(:,1)).*sin(pi.*p(:,2));

cval(:,1) =  pi.*cos(pi.*p(:,2)).*sin(pi.*p(:,1));
cval(:,2) = -pi.*cos(pi.*p(:,1)).*sin(pi.*p(:,2));