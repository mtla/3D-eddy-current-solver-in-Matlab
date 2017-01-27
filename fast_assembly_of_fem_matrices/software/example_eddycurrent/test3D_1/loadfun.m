function val = loadfun( p )

% LOADFUN
%   The load function of the PDE: curl curl u + u = f .
%
% SYNTAX:  val = loadfun( p )
%
% IN:   p      vector defining the points
%
% OUT:  val    values of the function f on points p
%

val = zeros(size(p));
val(:,3) = (2*pi^2 + 1).*sin(pi.*p(:,1)).*sin(pi.*p(:,2));