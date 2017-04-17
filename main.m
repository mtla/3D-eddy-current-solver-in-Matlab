DT = readMesh(strcat(pwd,'\meshes\example_mesh_3D.obj'));
tetramesh(DT); % plot mesh

buildStiffnessMatrix(DT)