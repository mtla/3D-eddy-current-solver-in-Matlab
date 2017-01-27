function favg = average( elems2nodes, B_K, b_K, B_K_det, f, qr)

% AVERAGE
%   Calculates the average of a function f (fun. handle).
%   Note that the user gives the quadrature rule to use for this
%   calculation!
%
% Here B_K, b_K, and B_K_det are given by the function
% 'affine_transformations(...)'.
%
% IN:   elems2nodes  elements by nodes
%       B_K          affine map matrices
%       b_K          affine map vectors
%       B_K_det      affine map determinants
%       f            function handle (i.e., the load function)
%       qr           the integration quadrature rule to use
%
% OUT:  favg         the average of f
%

if ( nargin==5 )
    qr=2;
end

dim    = size(B_K,1);
nelems = size(elems2nodes,1);

B_K_detA = abs(B_K_det);

% get the integration points and weighs on the ref element
[ip,w,nip] = intquad(qr,dim);

favg = zeros(nelems,1);

% calculate the integral of f elementwise
for i=1:nip
    %calculate F_K( integration points ) in the i'th int point
    F_K_ip = squeeze(amsv(B_K, ip(i,:)))' + b_K;
    
    %get the load function values at these transformed integration points
    fval = f(F_K_ip);
    
    favg = favg + w(i) .* B_K_detA .* fval;
end

% multiply integral values by areas of the elements to obtain the average
if ( dim==2 )
    favg = favg./(B_K_detA./2);
elseif ( dim==3 )
    favg = favg./(B_K_detA./6);
end
