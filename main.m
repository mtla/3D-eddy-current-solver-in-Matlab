msh = readMesh(strcat(pwd,'\meshes\example_mesh_3D.obj'));
dirichletNodes = readDirichletNodes(msh, strcat(pwd,'\meshes\example_mesh_3D_dirichlet.obj'));

% msh = readMesh(strcat(pwd,'\meshes\long_bar.obj'));
% dirichletNodes = readDirichletNodes(msh, strcat(pwd,'\meshes\long_bar_dirichlet.obj'));
figure(1)
tetramesh(msh); % plot mesh

permeability = 10;
permittivity = 10;

%reluctivity of each element [A/(Tm)]
reluctivity = 1/(pi*4e-7);
np = size(msh.Points,1);
%current density in each element [A/m^2]
% currentDensity = ones(np,1); % last element has a source current
% currentDensity(dirichletNodes) = 0;
%dirichletNodes;

[Sn, Se] = buildStiffnessMatrix(msh, permittivity);
[fn, fe] = buildLoadVector(msh);

freeNodes = setdiff(1:np, dirichletNodes); %nodes NOT in Dirichlet bnd

% calculating potentials in the free nodes
Afree = Sn(freeNodes,freeNodes) \ fn(freeNodes);
% NOTE: this is equivalent to Afree = inv(S) * f, but much faster

Aedges = Se \ fe;

%assembling solution in the entire region
A_total = zeros(np,1);
A_total(freeNodes) = Afree;

figure(2)
scatter3(msh.Points(:,1),msh.Points(:,2),msh.Points(:,3),A_total);
plot3(msh, Aedges);
% tetramesh(msh.ConnectivityList, msh.Points, Afree);
