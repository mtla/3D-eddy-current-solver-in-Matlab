verts = readObj(strcat(pwd,'\meshes\example_mesh_3D.obj'));
DT = delaunayTriangulation(verts);
tetramesh(DT); % plot mesh

tetrahedrons = DT.ConnectivityList;
tetrahedron = tetrahedrons(1,:);
points_list = DT.Points;

% test tetrahedron2matrix
tetrahedron2matrix(tetrahedron, points_list)