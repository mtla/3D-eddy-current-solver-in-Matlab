function [B,b] = get_ElementwiseMapping(msh, el)
%returns the mapping from the reference element to the global element el,
%such that g_global = B*x_ref + b
%
% Note that the mapping is now in matrix form, i.e.
% x = [(x_j-x_i) (x_k-x_i)]*[psi; eta] + x_i


B = [msh.p(:,msh.t(2,el))-msh.p(:,msh.t(1,el)) msh.p(:,msh.t(3,el))-msh.p(:,msh.t(1,el))];
b = msh.p(:,msh.t(1,el));
end