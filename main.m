DT = readMesh(strcat(pwd,'\meshes\example_mesh_3D.obj'));
tetramesh(DT); % plot mesh

tetrahedrons = DT.ConnectivityList;
tetrahedron = tetrahedrons(1,:);
points_list = DT.Points;

% test tetrahedron2matrix
tetrahedron2matrix(tetrahedron, points_list)