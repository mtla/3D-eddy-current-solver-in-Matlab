function [ sMatrix ] = buildStiffnessMatrix( DT )
% BUILDSTIFFNESMATRIX 
% This function inputs a delaunayTriangulation (struct), that is basically
% a mesh that has been divided into smaller tetrahedrons. It then
% calculates the stiffness matrix for this mesh and returns it.
%
% input:
%
%   delaunayTriangulation with properties:
%               Points: [n�3 double]
%     ConnectivityList: [m�4 double]
%          Constraints: [] (usually empty)
%
% output: [n�n double] matrix
%
    
    tetrahedrons = DT.ConnectivityList;
    vertices_list = DT.Points;
    sMatrix = zeros(max(max(tetrahedrons)));
    % get rid of the for loop. Matlab does not like them that much
    for row = 1:size(tetrahedrons, 1)
        tetrahedron = tetrahedrons(row, :);
        S = tetrahedron2Smatrix(tetrahedron, vertices_list);
        for i = 1:4
            for j = 1:4
                sMatrix(tetrahedron(i), tetrahedron(j)) = ...
                   sMatrix(tetrahedron(i), tetrahedron(j)) +  S(i,j);
            end
        end 
    end
end

function [ S_local ] = tetrahedron2Smatrix( tetrahedron, node_coordinates )
% TETRAHEDRON2MATRIX This functions takes an tetrahedron, along with a list
% of the coordinates of the nodes, as an input and outputs an 4x4 matrix
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
    % imagine it like a tetrahedron with the points:
    % P1 = [ 0 0 0 ]
    % P2 = [ 1 0 0 ]
    % P3 = [ 0 1 0 ]
    % P4 = [ 0 0 1 ]
    % If you want to know the volume of the tetrahedron and start
    % calculating from P1, the volume increases whether you go along the
    % x-, y- or z-axis. If you start from the point P2, the volume will
    % only increse if you move along the x-axis. P3, only y-axis etc.
    % This is also constant everywhere in the tetrahedron
    gradPhi_ref = [-1 -1 -1;1 0 0; 0 1 0;0 0 1]';
    w1 = 0.5; %integration weight for the single-point quadrature

    [B,~] = map2global(tetrahedron, node_coordinates);
    gradPhi = (B') \ gradPhi_ref; %gradients of shape functions of the GLOBAL element
    
    %assembling the element-contribution to the stiffness matrix
    %only upper triangular parts first
    S_local = zeros(4,4);
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



