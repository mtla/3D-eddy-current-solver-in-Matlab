function [ load_vector ] = buildLoadVector( DT )
% BUILDSTIFFNESMATRIX 
% This function inputs a delaunayTriangulation (struct), that is basically
% a mesh that has been divided into smaller tetrahedrons. It then
% calculates the stiffness matrix for this mesh and returns it.
%
% input:
%
%   delaunayTriangulation with properties:
%               Points: [n×3 double]
%     ConnectivityList: [m×4 double]
%          Constraints: [] (usually empty)
%
% output: [n×n double] matrix
%
    
    tetrahedrons = DT.ConnectivityList;
    vertices_list = DT.Points;
    load_vector = zeros(max(max(tetrahedrons)));
    % get rid of the for loop. Matlab does not like them that much
    for row = 1:size(tetrahedrons, 1)
        tetrahedron = tetrahedrons(row, :);
        S = tetrahedron2Lvector(tetrahedron, vertices_list);
        for i = 1:4
            for j = 1:4
                load_vector(tetrahedron(i), tetrahedron(j)) = ...
                   load_vector(tetrahedron(i), tetrahedron(j)) +  S(i,j);
            end
        end 
    end
end

function [ S_local ] = tetrahedron2Lvector( tetrahedron, node_coordinates )
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
%       load vector:
%           4x4 matrix with values according to the shape function
%
% Syntax: tetrahedron2Lvector(tetrahedron, nodes_coordinates)
%
    
    %reference shape functions expressed in polynomial basis
    Phi_ref = [1 -1 -1 -1;0 1 0 0;0 0 1 0;0 0 0 1]';
    
    xa = [0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105]; 
    ya = [0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685];
    za = [0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105];
    points = [xa;ya;za];
    w_quad = [.25 .25 .25]/6;
end



