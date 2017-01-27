function val = permeability( p )

% PERMEABILITY
%   The permeability mu of the PDE: curl mu^-1 curl u + kappa u = f .
%
% SYNTAX:  val = loadfun( p )
%
% IN:   p      vector defining the points
%
% OUT:  val    values of the function mu on points p
%

M   = size(p,1);
val = ones(M,1);
val(p(:,1)<0.5) = 10;