function [ sMatrixNodes, sMatrixEdges ] = buildStiffnessMatrix(msh, reluctivity)
% BUILDSTIFFNESMATRIX 
% This function inputs a delaunayTriangulation (struct), that is basically
% a mesh that has been divided into smaller tetrahedrons. It then
% calculates the stiffness matrix for this mesh and returns it.
%
% input:
%
%   delaunayTriangulation with properties:
%               Points: [n×3 double]
%     TetrahedronsByPoints: [m×4 double]
%          Constraints: [] (usually empty)
%
% output: [n×n double] matrix
%
    
    sMatrixNodes = zeros(msh.np());
    sMatrixEdges = zeros(msh.ne());
    % get rid of the for loop. Matlab does not like them that much
    % tID stands for tetrahedron ID. That is, where in the
    % msh.TetrahedronsByPoints array the tetrahedron is
    for tID = 1:msh.nt()
        % we can calculate the affine transformation 
        [B,~] = map2global(msh, tID);
        S = points2Smatrix(B);
        tetrahedron = msh.tetrahedron2points(tID)
        for i = 1:4
            for j = 1:4
                sMatrixNodes(tetrahedron(i), tetrahedron(j)) = ...
                   sMatrixNodes(tetrahedron(i), tetrahedron(j)) +  S(i,j);
            end
        end
        edges = msh.tetrahedron2edges(tID);
        Sedges = edges2Smatrix(B);
    end
    sMatrixNodes = sMatrixNodes * reluctivity;
%     for row = 1:size(msh.nt(),1)
%         S = edges2Smatrix(msh, row);
%         tetrahedron = msh.TetrahedronsByEdges(row,:);
%         for i = 1:4
%             for j = 1:4
%                 sMatrixEdges(tetrahedron(i), tetrahedron(j)) = ...
%                    sMatrixEdges(tetrahedron(i), tetrahedron(j)) +  S(i,j);
%             end
%         end
%     end
    
end

function [ S_local ] = points2Smatrix(B)
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
    % the gradient of a 3D tetrahedron
    % imagine it like a tetrahedron with the points and values:
    % P1 = 1 - x - y - z
    % P2 = x
    % P3 = y
    % P4 = z
    % When we take the gradient of this, we get
    gradPhi_ref = [-1 -1 -1;1 0 0; 0 1 0;0 0 1]';
    w1 = 0.5; %integration weight for the single-point quadrature

%     [B,~] = map2global(msh, node_coordinates);
    gradPhi = (B') \ gradPhi_ref %gradients of shape functions of the GLOBAL element
    
    %assembling the element-contribution to the stiffness matrix
    %only upper triangular parts first
    S_local = zeros(4);
    S_local(1,2) = gradPhi(:,1)' * gradPhi(:,2);
    S_local(1,3) = gradPhi(:,1)' * gradPhi(:,3);
    S_local(1,4) = gradPhi(:,1)' * gradPhi(:,4);
    S_local(2,3) = gradPhi(:,2)' * gradPhi(:,3);
    S_local(2,4) = gradPhi(:,2)' * gradPhi(:,4);
    S_local(3,4) = gradPhi(:,3)' * gradPhi(:,4);
    S_local = S_local + S_local'; % S_local is symmetrical so we can get the lower part by summing its transpose
    
    % calculate the diagonal
    S_local(1,1) = gradPhi(:,1)' * gradPhi(:,1);
    S_local(2,2) = gradPhi(:,2)' * gradPhi(:,2);
    S_local(3,3) = gradPhi(:,3)' * gradPhi(:,3);
    S_local(4,4) = gradPhi(:,4)' * gradPhi(:,4);
    % Calculates the actual contribution of the tetrahedron by using single
    % point quadrature
    % TODO: add the contribution of the reluctivity
    S_local = w1 * S_local * abs(det(B));
end

function [ S_local ] = edges2Smatrix(B)
%     edges
    S_local = zeros(6);
%     [B,~] = map2global(msh, edges);
end
