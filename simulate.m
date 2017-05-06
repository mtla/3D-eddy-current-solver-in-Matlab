[filename, filefolder] = uigetfile({'*.obj','Object Files';'*.txt','Raw Text File'}, 'Read obj-file');
mesh_path = strcat(filefolder,filename);

[filename, filefolder] = uigetfile({'*.obj','Object Files';'*.txt','Raw Text File'}, 'Read obj-file');
dirichlet_path = strcat(filefolder,filename);

msh = readMesh(mesh_path);
dirichletNodes = readDirichletNodes(msh, dirichlet_path);

prompt = {'Permittivity in element:','Permeability in element:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {num2str(permittivity), num2str(permeability)};
constants = inputdlg(prompt,dlg_title,num_lines,defaultans);
[permittivity, permeability] = constants{:};

permittivity = num2str(permittivity);
% plot mesh
figure(1)
tetramesh(msh);

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
scatter3(msh.Points(:,1),msh.Points(:,2),msh.Points(:,3),abs(A_total)*10^9+1);
plot3(msh, Aedges);

[filename, filefolder] = uisavefile({'*.csv','Comma separated file';'*.txt','Raw Text File';'*.vtk','Binary Paraview file'}, 'Save output');
output_path = strcat(filefolder,filename);