function [ LvectorNodes, LvectorEdges ] = buildLoadVector(msh) 
% This function inputs a delaunayTriangulation (struct), that is basically
% a mesh that has been divided into smaller tetrahedrons. It then
% calculates the load vector for this mesh and returns it.
%
% input:
%
%   msh with properties:
%                   Points: [n×3 double]
%     TetrahedronsByPoints: [m×4 double]
%              Constraints: [] (usually empty)
%
% output: [n×1 double] matrix
%

    LvectorNodes = zeros(msh.np(), 1);
    LvectorEdges = zeros(msh.ne(), 1);
    
    %reference shape functions expressed in polynomial basis
    Phi_ref = [1 -1 -1 -1;0 1 0 0;0 0 1 0;0 0 0 1]';
    
    % reference shape function for edges expressed in polynomial basis
    ePhi_ref = [ Phi_ref(:,2)-Phi_ref(:,1) ...
                 Phi_ref(:,3)-Phi_ref(:,1) ...
                 Phi_ref(:,4)-Phi_ref(:,1) ...
                 Phi_ref(:,3)-Phi_ref(:,2) ...
                 Phi_ref(:,4)-Phi_ref(:,3) ...
                 Phi_ref(:,2)-Phi_ref(:,4) ];
                 
        
    
    % get integration points for reference tetrahedron
    [quad_points, w_quad] = inttet(2);
    % Transpose matrixes because inttet thinks in a translated way
    quad_points = quad_points';
    w_quad = w_quad';
    
%     [values, ~] = basis_Nedelec0(quad_points');
    
    for row = 1:msh.nt()
        tetrahedron = msh.TetrahedronsByPoints(row, :);
        edges = msh.tetrahedron2edges(row);
%         L = tetrahedron2Lvector(msh);
        [B,~] = map2global(msh, row);

        %looping over the integration points
        for k_quad = 1:4
            w1 = w_quad(k_quad);
            psi = [1 quad_points(1,k_quad) quad_points(2,k_quad) quad_points(3,k_quad)]; %polynomial basis at the integration point
            Phi = psi * Phi_ref; %shape function values at the integration point
%             val = w1 * Phi' * abs(det(B));
            Phi_edges = psi * ePhi_ref;
            B_det = abs(det(B));
            LvectorNodes(tetrahedron') = LvectorNodes(tetrahedron') + ...
                w1 * Phi' * B_det;
            LvectorEdges(abs(edges')) = LvectorEdges(abs(edges')) + ...
                w1 * Phi_edges' * B_det; 
        end
%         for k_quad = 1:6
    end
end