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
val(p(:,2)<0.5,1) = 1;