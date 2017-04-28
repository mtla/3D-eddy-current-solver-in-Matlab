DT = readMesh(strcat(pwd,'\meshes\long_bar.obj'));
dirichletNodes = readDirichletNodes(strcat(pwd,'\meshes\example_mesh_3D_dirichlet.obj'), DT);
tetramesh(DT); % plot mesh

%reluctivity of each element [A/(Tm)]
reluctivity = (pi*4e-7);
np = size(DT.Points,1);
%current density in each element [A/m^2]
currentDensity = ones(np,1); % last element has a source current [ 0 1 0 1 0 1]
currentDensity(dirichletNodes) = 0;
%dirichletNodes;

S = buildStiffnessMatrix(DT, reluctivity);
f = buildLoadVector(DT, currentDensity);


freeNodes = setdiff(1:np, DirichletNodes); %nodes NOT in Dirichlet bnd

%calculating potentials in the free nodes
Afree = S(freeNodes,freeNodes) \ f(freeNodes);
% NOTE: this is equivalent to Afree = inv(S) * f, but much faster

%assembling solution in the entire region
A_total = zeros(Np,1);
A_total(freeNodes) = Afree;

figure
scatter3(DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),abs(A_total)/100+1);
% tetramesh(DT.ConnectivityList, DT.Points, Afree);
