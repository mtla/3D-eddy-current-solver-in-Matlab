function [ sMatrix ] = buildStiffnessMatrix( tetrahedrons, vertices_list )
%BUILDSTIFFNESMATRIX Summary of this function goes here
%   Detailed explanation goes here

    sMatrix = zeros(max(max(tetrahedrons)));
    for row = 1:size(tetrahedrons, 1)
        tetrahedron = tetrahedrons(row, :);
        S = tetrahedron2matrix(tetrahedron, vertices_list);
        for i = 1:4
            for j = 1:4
                sMatrix(tetrahedron(i), tetrahedron(j)) = ...
                   sMatrix(tetrahedron(i), tetrahedron(j)) +  S(i,j);
            end
        end 
    end
end

