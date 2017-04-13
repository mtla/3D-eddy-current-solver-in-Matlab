readObj('D:\Henrik\Programs\eddy-currents-fem\meshes\example_mesh_3D.obj');
verts = ans;
DT = delaunayTriangulation(verts.vertices);
tetramesh(DT); % plot mesh

tetrahedrons = DT.ConnectivityList;
tetrahedron = tetrahedrons(1,:);
points_list = DT.Points;

% test tetrahedron2matrix
tetrahedron2matrix(tetrahedron, points_list)