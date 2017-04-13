function [ tetrahedron_matrix ] = tetrahedron2matrix( tetrahedron, node_coordinates )
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
    
    tetrahedron_matrix = zeros(4,4);
%     volume = det([ones(4,1) x' y' z'])/6;
    [B,~] = map2global(tetrahedron, node_coordinates)
    gradPhi = (B') \ gradPhi_ref %gradients of shape functions of the GLOBAL element
    
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
    S_local
    
%     for i = 1:4
%         for j = 1:4
%             % Determines the value of the specific S_ij
%             % with the help of the shapefunction
%             % 
%             % CURRENT SHAPE FUNCTION IS NOT THE CORRECT ONE
% %             tetrahedron_matrix(i,j) = (x(i)-x(j))+(y(j)-y(i));
%         end
%     end
end

function [ B, b ] = map2global( tetrahedron, node_coordinates ) 
% Maps a single tetrahedron to the global unit tetrahedron so that
% g_global = B*x_ref + b

    x = zeros(1,4);
    y = zeros(1,4);
    z = zeros(1,4);
    for n = 1:4
        x(n) = node_coordinates(tetrahedron(n),1);
        y(n) = node_coordinates(tetrahedron(n),2);
        z(n) = node_coordinates(tetrahedron(n),3);
    end
    m_tetrahedron = [x;y;z];
    B = [m_tetrahedron(:,4)-m_tetrahedron(:,1) m_tetrahedron(:,3)-m_tetrahedron(:,1) m_tetrahedron(:,2)-m_tetrahedron(:,1)];
    b = m_tetrahedron(:,1);
end

