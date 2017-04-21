DT = readMesh(strcat(pwd,'\meshes\example_mesh_3D.obj'));
dirichletNodes = readDirichletNodes(strcat(pwd,'\meshes\example_mesh_3D_dirichlet.obj'), DT);
tetramesh(DT); % plot mesh

% buildStiffnessMatrix(DT)