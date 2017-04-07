function [S,f] = calculate_Arrays(msh, reluctivityInElement,sourceInElement)
%calculate_Arrays returns the stiffness matrix and load vector
% [S,f] = calculate_Arrays(msh, reluctivityInElement,sourceInElement) returns
% the stiffness matrix S and the load vector f for a magnetostatic problem.
% NOTE: boundary conditions are ignored here.
%   
%   Input arguments:
%   - msh: a struct containing the following fields
%       .p: a 2xNp matrix, each column of which contains the coordinates of
%       a node
%       .t: a 3xNp matrix, each column of which defines a triangle
%   - reluctivityInElement: a vector of Ne entries, containing the material
%   reluctivity inside each element (assumed constant inside any single
%   element)
%   - sourceInElement: a vector of Ne entries, containing the current
%   density inside each element

S = calculate_StiffnessMatrix(msh, reluctivityInElement);
f = calculate_LoadVector(msh, sourceInElement);

end

function S = calculate_StiffnessMatrix(msh, reluctivityInElement)
%returns the stiffness matrix with the entries
%S_{ij} = \Int{ reluctivity(x) * (\grad Phi_i) \cdot (\grad Phi_j)

%The integrals are calculated with a
%single-point Gaussian quadrature 
%(sufficient for first-order shape functions and mesh)

Ne = size(msh.t, 2); %number of elements in the mesh
Np = size(msh.p, 2); %number of nodes

S = zeros(Np, Np); %initializing the matrix
%NOTE: normally, a sparse matrix would be assembled. However, for this kind
%of a small demonstrative problem, dealing with standard (dense) matrices is both
%faster and easier to follow

%gradients of the reference element shape functions (constant over the element)
%(also hard-coded for extra ugliness)
gradPhi_ref = [-1 -1;1 0; 0 1]';

w1 = 0.5; %integration weight for the single-point quadrature

%looping over elements
for k_element = 1:Ne
    [B,~] = get_ElementwiseMapping(msh, k_element);    
    gradPhi = (B') \ gradPhi_ref; %gradients of shape functions of the GLOBAL element
    
    %assembling the element-contribution to the stiffness matrix
    %only upper triangular parts first
    S_local = zeros(3,3);
    S_local(1,1) = gradPhi(:,1)' * gradPhi(:,1);
    S_local(2,2) = gradPhi(:,2)' * gradPhi(:,2);
    S_local(3,3) = gradPhi(:,3)' * gradPhi(:,3);
    S_local = S_local / 2; %why? --> see row 61
    
    S_local(1,2) = gradPhi(:,1)' * gradPhi(:,2);
    S_local(1,3) = gradPhi(:,1)' * gradPhi(:,3);
    S_local(2,3) = gradPhi(:,2)' * gradPhi(:,3);
    S_local = S_local + S_local'; %getting the lower diagonal parts
    
    %adding the contribution of k_element to the global matrix
    indices = msh.t(:,k_element); %nodes making this element
    
    %calculating the contribution with a single-point quadrature
    S(indices, indices) = S(indices, indices) + ...
        w1 * reluctivityInElement(k_element) * S_local * abs(det(B));
end

end

function f = calculate_LoadVector(msh, sourceInElement)
%returns the load vector
% 
% 3-point quadrature is used to compute the integral.

Ne = size(msh.t, 2);
Np = size(msh.p, 2);

f = zeros(Np, 1);

%reference shape functions expressed in polynomial basis
Phi_ref = [1 -1 -1;0 1 0;0 0 1]';

%values of reference shape functions in the quadrature point [1/3 1/3];
x_quad = [1/6 2/3 1/6; 1/6 1/6 2/3];
w_quad = [1/6 1/6 1/6];

for k_element = 1:Ne
    [B,~] = get_ElementwiseMapping(msh, k_element);
    
    indices = msh.t(:,k_element);
	
	%looping over the integration points
	for k_quad = 1:3
		w1 = w_quad(k_quad);
		psi = [1 x_quad(1,k_quad) x_quad(2,k_quad)]; %polynomial basis at the integration point
		Phi = psi * Phi_ref; %shape function values at the integration point
		
		f(indices) = f(indices) + ...
			w1 * sourceInElement(k_element) * Phi' * abs(det(B));
	end
end

end