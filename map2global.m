function [ B, b ] = map2global(DT, element) 
% Maps a single tetrahedron to the global unit tetrahedron so that
% g_global = B*x_ref + b
    coordinates = DT.Points;
    t = DT.ConnectivityList(element,:); % tetrahedron
    % build the tetrahedron matrix so that
    % matrix = [x_1...x_4 ; y_1...y_4 ; z_1...z_4]
    gm = zeros(3,4); % global matrix
    for n = 1:4
        gm(:,n) = coordinates(t(n),:)';
    end
    B = [gm(:,2)-gm(:,1) gm(:,3)-gm(:,1) gm(:,4)-gm(:,1)];
    b = gm(:,1);
end