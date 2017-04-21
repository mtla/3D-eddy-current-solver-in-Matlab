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
        S = tetrahedron2matrix(tetrahedron, vertices_list);
        for i = 1:4
            for j = 1:4
                load_vector(tetrahedron(i), tetrahedron(j)) = ...
                   load_vector(tetrahedron(i), tetrahedron(j)) +  S(i,j);
            end
        end 
    end
end

