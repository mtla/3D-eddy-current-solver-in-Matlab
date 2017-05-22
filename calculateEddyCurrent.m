function [ J ] = calculateEddyCurrent( msh, T )
%CALCULATEEDDYCURRENT Summary of this function goes here
%   Detailed explanation goes here
    J = zeros(size(msh.nt(),1),3);
    
    for i = 1:msh.nt()
        tetrahedron = msh.TetrahedronsByEdges(i, :)
        
        [B,~] = map2global(msh, i);
        [integration_points, weights] = inttet(1);
        [~, curl_ref] = basis_Nedelec0(integration_points);

        curlPhi = B * curl_ref / det(B);
        J(i,:) = T(abs(tetrahedron(1))) * curlPhi(:,1)' + ...
            T(abs(tetrahedron(2))) * curlPhi(:,2)' + ...
            T(abs(tetrahedron(3))) * curlPhi(:,3)' + ...
            T(abs(tetrahedron(4))) * curlPhi(:,4)' + ...
            T(abs(tetrahedron(5))) * curlPhi(:,5)' + ...
            T(abs(tetrahedron(6))) * curlPhi(:,6)';
            
    end

end

