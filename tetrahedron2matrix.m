function [ tetrahedron_matrix ] = tetrahedron2matrix( tetrahedron, nodes_coordinates )
% TETRAHEDRON2MATRIX This functions takes an tetrahedron, along with a list
% of the coordinates of the nodes, as an input and outputs an 4x4 matrix
%
% Input:
%       tetrahedron: 
%           array with for elements
%       nodes_coordinates: 
%           3xn matrix with coordinates for each vertice in
%           the mesh
%
% Output:
%       tetrahedron_matrix:
%           4x4 matrix with values according to the shape function
%
% Syntax: tetrahedron2matrix(tetrahedron, nodes_coordinates)
%

    tetrahedron_matrix = zeros(4,4);
    x = zeros(1,4);
    y = zeros(1,4);
    z = zeros(1,4);
    for n = 1:4
        x(n) = nodes_coordinates(tetrahedron(n),1);
        y(n) = nodes_coordinates(tetrahedron(n),2);
        z(n) = nodes_coordinates(tetrahedron(n),3);
    end
    for i = 1:4
        for j = 1:4
            % Determines the value of the specific S_ij
            % with the help of the shapefunction
            % 
            % CURRENT SHAPE FUNCTION IS NOT THE CORRECT ONE
            if (i ~= j)
                tetrahedron_matrix(i,j) = (x(i)-x(j))+(y(j)-y(i));
            end
        end
    end
end

