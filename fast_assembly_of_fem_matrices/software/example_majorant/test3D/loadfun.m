function val = loadfun( p )

% LOADFUN
%   The load function of the PDE: - div grad u = f .
%
% SYNTAX:  val = loadfun( p )
%
% IN:   p      vector defining the points
%
% OUT:  val    values of the function f on points p
%

val = - 2.*p(:,1).*p(:,2).*(p(:,1) - 1).*(p(:,2) - 1) - 2.*p(:,1).*p(:,3).*(p(:,1) - 1).*(p(:,3) - 1) - 2.*p(:,2).*p(:,3).*(p(:,2) - 1).*(p(:,3) - 1);