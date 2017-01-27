function [residual,dual] = majorant( elems2nodes, elems, B_K, b_K, B_K_det, signs, f, u_h, y_h )

% MAJORANT
%   Calculates the majorant's residual and dual values elementwise with
%   the given approximations u_h and y_h.
%
% SYNTAX:   maj = majorant( elements, elems2faces, B_K, b_K, B_K_det, signs, f, u_h, y_h )
%
% IN:   elems2nodes  elements by nodes
%       elems        elements by: edges in 2D / faces in 3D
%       B_K          affine map matrices
%       b_K          affine map vectors
%       B_K_det      affine map determinants
%       signs        RT basis function signs per element
%       u_h          potential approximation coefficients
%       y_h          flux approximation coefficients
%
% OUT:  residual     the value of the majorant's residual term (elementwise, squared)
%       dual         the value of the majorant's dual term (elementwise, squared)
%

dim      = size(B_K,1);
nelems   = size(elems2nodes,1);
B_K_detA = abs(B_K_det);

residual = zeros(nelems,1);
dual     = zeros(nelems,1);

% calculate the residual part
%----------------------------
qrR = 6;                           % which integration quadrature to use
[ipR,wR,nipR] = intquad(qrR,dim);  % R as in residual

[~,y_h_dval] = solution_RT0(elems,B_K,B_K_det,signs,y_h,1);  % flux approx. div (constant per element)

for i=1:nipR
    %calculate F_K( integration points ) in the i'th int point
    F_K_ip = squeeze(amsv(B_K, ipR(i,:)))' + b_K;
    
    %get the exact solution values at the transformed integration points
    fval = f(F_K_ip);
    
    %the residual
    R = fval + y_h_dval;
    
    residual = residual + wR(i) * B_K_detA .* R.^2;
end

% calculate the dual part
%------------------------
qrD = 2;                           % which integration quadrature to use
[~,wD,nipD] = intquad(qrD,dim);    % D as in dual

[~      ,u_h_gval] = solution_P1(elems2nodes,B_K,u_h,qrD);           % potential approx. grad
[y_h_val,~       ] = solution_RT0(elems,B_K,B_K_det,signs,y_h,qrD);  % flux approx.

for i=1:nipD
    %the dual
    D = y_h_val(:,:,i) - u_h_gval(:,:,i);
    
    dual = dual + wD(i) * B_K_detA .* sum( D.^2, 2 );
end