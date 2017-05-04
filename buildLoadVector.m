function [ load_vector ] = buildLoadVector(msh, currentDensity) 
% This function inputs a delaunayTriangulation (struct), that is basically
% a mesh that has been divided into smaller tetrahedrons. It then
% calculates the load vector for this mesh and returns it.
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
    np = size(msh.Points, 1); % number of points
    ne = size(msh.TetrahedronsByPoints, 1); % number of tetrahedrons
    load_vector = zeros(np,1);
    % get rid of the for loop. Matlab does not like them that much
    
    %reference shape functions expressed in polynomial basis
    Phi_ref = [1 -1 -1 -1;0 1 0 0;0 0 1 0;0 0 0 1]';

    xa = [0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105]; 
    ya = [0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685];
    za = [0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105];
    quad_points = [xa;ya;za];
    w_quad = [.25 .25 .25 .25]/6;
    
    for row = 1:ne
        tetrahedron = msh.TetrahedronsByPoints(row, :);
        [B,~] = map2global(msh, row);
%         L = tetrahedron2Lvector(msh);

        %looping over the integration points
        for k_quad = 1:4
            w1 = w_quad(k_quad);
            psi = [1 quad_points(1,k_quad) quad_points(2,k_quad) quad_points(3,k_quad)]; %polynomial basis at the integration point
            Phi = psi * Phi_ref; %shape function values at the integration point
%             val = w1 * Phi' * abs(det(B));
            load_vector(tetrahedron') = load_vector(tetrahedron') + ...
                w1 * Phi' * abs(det(B)) * currentDensity(row);
        end
    end
end