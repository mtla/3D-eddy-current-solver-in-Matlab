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
disp('Calculating Stiffness matrix...')
[Sn, Se] = buildStiffnessMatrix(msh, permittivity);
disp('Done!')
disp('Calculating Load vector...')
[fn, fe] = buildLoadVector(msh);
disp('Done!')
disp('Calculating C matrix...')
C = buildCMatrix(msh, permeability);
disp('Done!')

% freeNodes = setdiff(1:np, dirichletNodes); %nodes NOT in Dirichlet bnd

% [T, Omega] = [Se C; C' Sn] \ [fe; fn]; 
% [T, Omega] = inv([Se C'; C Sn]) * [fe; fn]

T = Se \ fe + C \ fn;
Omega = C' \ fe + Sn \ fn;

J = calculateEddyCurrent(msh,T);
centroid = msh.getCentroid();

% % calculating potentials in the free nodes
% Afree = Sn(freeNodes,freeNodes) \ fn(freeNodes);
% % NOTE: this is equivalent to Afree = inv(S) * f, but much faster
% 
% T = Se \ fe;
% 
% %assembling solution in the entire region
% Omega = zeros(np,1);
% Omega(freeNodes) = Afree;
% 
% msh.setPointValues(Omega);
% msh.setEdgeValues(T);
% 
figure(2)
% % scatter3(msh.Points(:,1),msh.Points(:,2),msh.Points(:,3),A_total);
scatter3(msh);
plot3(msh, T);
% tetramesh(msh.ConnectivityList, msh.Points, Afree);
quiver3(centroid(:,1),centroid(:,2),centroid(:,3), ...
    J(:,1),J(:,2),J(:,3))

% writeResults(msh)

disp('Done!')