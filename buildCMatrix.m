function [ C ] = buildCMatrix(msh, permeability)
% BUILDSTIFFNESMATRIX 
% This function inputs a delaunayTriangulation (struct), that is basically
% a mesh that has been divided into smaller tetrahedrons. It then
% calculates the stiffness matrix for this mesh and returns it.
%
% input:
%
%   mesh with properties:
%               Points: [n×3 double]
%     TetrahedronsByPoints: [m×4 double]
%          Constraints: [] (usually empty)
%
% output: [n×m double] matrix
%
    
    C = zeros(msh.np(),msh.ne());
    % get rid of the for loop. Matlab does not like them that much
    % tID stands for tetrahedron ID. That is, where in the
    % msh.TetrahedronsByPoints array the tetrahedron is
    for tID = 1:msh.nt()
        % we can calculate the affine transformation 
        [B,~] = map2global(msh, tID);
        C_local = B2Cmatrix(B) * permeability;
        tetrahedron = msh.tetrahedron2points(tID);
        edges = msh.tetrahedron2edges(tID);
        for i = 1:6
            for j = 1:4
                C(edges(i), tetrahedron(j)) = ...
                   C(edges(i), tetrahedron(j)) +  C_local(i,j);
            end
        end
    end
    C = C * permeability;
end

function [ C_local ] = B2Cmatrix(B)
% TETRAHEDRON2MATRIX This functions takes an DelaunayTriangulation (struct)
% along with the position of the tetrahedron (in the struct)
%
% Input:
%       tetrahedron: 
%           array with four elements
%       nodes_coordinates: 
%           3xn matrix with coordinates for each vertice in
%           the mesh
%
% Output:
%       Local S-matrix:
%           4x4 matrix with values according to the shape function
%
% Syntax: tetrahedron2Smatrix(tetrahedron, nodes_coordinates)
%

    [integration_points, weights] = inttet(1);
    [f_ref, ~] = basis_Nedelec0(integration_points);
    
    W = (B')^-1 * f_ref;
    % imagine it like a tetrahedron with the points and values:
    % P1 = 1 - x - y - z
    % P2 = x
    % P3 = y
    % P4 = z
    % When we take the gradient of this, we get
    gradPhi_ref = [-1 -1 -1;1 0 0; 0 1 0;0 0 1]';
    w1 = 0.5; %integration weight for the single-point quadrature
    
    % gradients of shape functions of the GLOBAL element
    gradPhi = (B') \ gradPhi_ref;
    % this is actually equal to (B')^-1 * gradPhi_ref
    % but is much faster
    
    %assembling the element-contribution to the stiffness matrix
    %only upper triangular parts first
    C_local = zeros(6,4);
    for i = 1:6
        C_local(i,1) = W(:,i)' * gradPhi(:,1);
        C_local(i,2) = W(:,i)' * gradPhi(:,2);
        C_local(i,3) = W(:,i)' * gradPhi(:,3);
        C_local(i,4) = W(:,i)' * gradPhi(:,4);
    end
    
    C_local = w1 * C_local;
end
