function [ S_local ] = tetrahedron2matrix( tetrahedron, node_coordinates )
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
%       tetrahedron_matrix:
%           4x4 matrix with values according to the shape function
%
% Syntax: tetrahedron2matrix(tetrahedron, nodes_coordinates)
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

function [ B, b ] = map2global( tetrahedron, node_coordinates ) 
% Maps a single tetrahedron to the global unit tetrahedron so that
% g_global = B*x_ref + b

    % build the tetrahedron matrix so that
    % matrix = [x_1...x_4 ; y_1...y_4 ; z_1...z_4]
    m_tetrahedron = zeros(3,4);
    for n = 1:4
        m_tetrahedron(1,n) = node_coordinates(tetrahedron(n),1);
        m_tetrahedron(2,n) = node_coordinates(tetrahedron(n),2);
        m_tetrahedron(3,n) = node_coordinates(tetrahedron(n),3);
    end
    B = [m_tetrahedron(:,2)-m_tetrahedron(:,1) m_tetrahedron(:,3)-m_tetrahedron(:,1) m_tetrahedron(:,4)-m_tetrahedron(:,1)];
    b = m_tetrahedron(:,1);
end

