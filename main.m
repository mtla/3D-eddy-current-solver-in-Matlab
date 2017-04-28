DT = readMesh(strcat(pwd,'\meshes\example_mesh_3D.obj'));
dirichletNodes = readDirichletNodes(strcat(pwd,'\meshes\example_mesh_3D_dirichlet.obj'), DT)
% tetramesh(DT); % plot mesh

%reluctivity of each element [A/(Tm)]
reluctivity = (pi*4e-7);

%current density in each element [A/m^2]
currentDensity = [ zeros(size(DT.ConnectivityList,1) - 1,1)' 1]; % last element has a source current
%dirichletNodes;

S = buildStiffnessMatrix(DT, reluctivity);
f = buildLoadVector(DT, currentDensity);

freeNodes = setdiff(1:size(DT.Points,1), DirichletNodes); %nodes NOT in Dirichlet bnd

%calculating potentials in the free nodes
Afree = S(freeNodes,freeNodes) \ f(freeNodes);
% NOTE: this is equivalent to Afree = inv(S) * f, but much faster
