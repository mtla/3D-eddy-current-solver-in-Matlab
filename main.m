read_obj
verts = ans;
DT = delaunayTriangulation(verts.vertices);
tetramesh(DT);

tetrahedrons = DT.ConnectivityList;
tetrahedron = tetrahedrons(1,:);
points_list = DT.Points;

% test tetrahedron2matrix
tetrahedron2matrix(tetrahedron, points_list)