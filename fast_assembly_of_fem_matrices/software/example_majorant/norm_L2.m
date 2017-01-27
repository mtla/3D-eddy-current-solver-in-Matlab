function fnorm = norm_L2( elems2nodes, B_K, b_K, B_K_det, f )

% NORM_L2
%   Calculates the L2-norm of a function f (fun. handle)
%
% Here B_K, b_K, and B_K_det are given by the function
% 'affine_transformations(...)'.
%
% IN:   elems2nodes  elements by nodes
%       B_K          affine map matrices
%       b_K          affine map vectors
%       B_K_det      affine map determinants
%       f            function handle (i.e., the load function)
%
% OUT:  norm         the L2-norm of f (squared)
%

dim    = size(B_K,1);
nelems = size(elems2nodes,1);

B_K_detA = abs(B_K_det);

% get the integration points and weighs on the ref element
[ip,w,nip] = intquad(6,dim);

fnorm = zeros(nelems,1);

for i=1:nip
    %calculate F_K( integration points ) in the i'th int point
    F_K_ip = squeeze(amsv(B_K, ip(i,:)))' + b_K;
    
    %get the load function values at these transformed integration points
    fval = f(F_K_ip);
    
    fnorm = fnorm + w(i) .* B_K_detA .* fval.^2;
end

fnorm = sum(fnorm);