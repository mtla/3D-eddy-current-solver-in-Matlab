function bfaces2nodes = get_boundary_faces(elems2faces,faces2nodes)

% GETBOUNDARYFACES
%   Calculates the boundary faces by their nodes for 3D mesh. In 2D this
%   data is calculated by 'refinement_uniform()'.
%
% IN:  elems2faces      elements by their faces
%      face2nodes       faces by their nodes
%
% OUT: bfaces2nodes     boundary faces by their nodes
%

[A1,I1] = sort(elems2faces(:,1));
[A2,I2] = sort(elems2faces(:,2));
[A3,I3] = sort(elems2faces(:,3));
[A4,I4] = sort(elems2faces(:,4));

nfaces = max(max(elems2faces));
E = zeros(nfaces,8);

E(A1,1) = I1;
E(A2,2) = I2;
E(A3,3) = I3;
E(A4,4) = I4;

% If the same face is listed in the same row of 'elems2faces' more than,
% once it will simply be missed! Because of this we have to insert the
% following dummy variables in order to determine the boundary faces.
ind1 = (diff(A1) == 0);
ind2 = (diff(A2) == 0);
ind3 = (diff(A3) == 0);
ind4 = (diff(A4) == 0);

E(A1(ind1),5) = 1;
E(A2(ind2),6) = 1;
E(A3(ind3),7) = 1;
E(A4(ind4),8) = 1;

% final sorting
E = sort(E,2,'descend');

% Get boundary nodes by first examining which columns in E
% have only one nonzero element, meaning that this face is
% related to only one single tetra, which means it is on the
% boundary of the domain. Since faces are defined by their nodes,
% we have the boundary nodes too.
ind = (E(:,2) == 0);
bfaces2nodes = faces2nodes(ind,:);

end