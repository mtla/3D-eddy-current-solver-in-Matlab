function [ sMatrix ] = buildStiffnessMatrix( tetrahedrons, vertices_list )
%BUILDSTIFFNESMATRIX Summary of this function goes here
%   Detailed explanation goes here

    sMatrix = zeros(size(tetrahedrons, 1));
%     for row = 1:size(tetrahedrons, 1)
%         tetrahedron = tetrahedrons(row, :);
%         S = tetrahedron2matrix(tetrahedron, vertices_list);
%     end
%     n = tetrahedrons(1,:)(end);
%     sMatrix = [n n];
end

