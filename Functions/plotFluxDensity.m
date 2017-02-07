function [] = plotFluxDensity(a, msh)
%plotFluxDensity plots the flux density amplitude as color.
% 
% plotFluxDensity(a, msh) plots the flux density amplitude, determined by
% the vector potential a and the mesh msh.

%getting the flux density in each element
Babs = getFluxDensity(a, msh);
Babs = repmat(Babs', 3, 1); %legacy fix for Matlab <2015

%describing the elements as polygons for Matlab
X = zeros(3, size(msh.t,2));
Y = X;

for kn = 1:3
    X(kn,:) = msh.p(1, msh.t(kn,:));
    Y(kn,:) = msh.p(2, msh.t(kn,:));
end

%plotting
fill(X,Y, Babs, 'Linestyle', 'none');
colormap('jet');
colorbar;
axis tight;

end

function Babs = getFluxDensity(a, msh)
Ne = size(msh.t, 2); %number of elements

Babs = zeros(Ne, 1);

gradPhi_ref = [-1 -1;1 0; 0 1]'; %ref. shapefun gradients again
for ke = 1:Ne
    [B,~] = get_ElementwiseMapping(msh, ke);
    indices = msh.t(:, ke);
    
    gradA = sum( bsxfun(@times, (B') \ gradPhi_ref, a(indices)'), 2);
    Babs(ke) = norm(gradA);    
end

end