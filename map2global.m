function [ B, b ] = map2global( tetrahedron, node_coordinates ) 
% Maps a single tetrahedron to the global unit tetrahedron so that
% g_global = B*x_ref + b

    % build the tetrahedron matrix so that
    % matrix = [x_1...x_4 ; y_1...y_4 ; z_1...z_4]
    m_tetrahedron = zeros(3,4);
    for n = 1:4
        m_tetrahedron(:,n) = node_coordinates(tetrahedron(n),:)';
    end
    B = [m_tetrahedron(:,2)-m_tetrahedron(:,1) m_tetrahedron(:,3)-m_tetrahedron(:,1) m_tetrahedron(:,4)-m_tetrahedron(:,1)];
    b = m_tetrahedron(:,1);
end