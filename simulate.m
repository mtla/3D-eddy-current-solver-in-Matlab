[filename, filefolder] = uigetfile({'*.obj','Object Files';'*.txt','Raw Text File'}, 'Read obj-file');
mesh_path = strcat(filefolder,filename);

[filename, filefolder] = uigetfile({'*.obj','Object Files';'*.txt','Raw Text File'}, 'Read obj-file');
dirichlet_path = strcat(filefolder,filename);

msh = readMesh(mesh_path);
dirichletNodes = readDirichletNodes(msh, dirichlet_path);

if (exist('permittivity','var')==0)
    permittivity = 10;
end
if (exist('permeability','var')==0)
    permeability = 10;
end

prompt = {'Permittivity in element:','Permeability in element:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {num2str(permittivity), num2str(permeability)};
constants = inputdlg(prompt,dlg_title,num_lines,defaultans);
[permittivity, permeability] = constants{:};

permittivity = str2double(permittivity);
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
scatter3(msh, A_total);
plot3(msh, Aedges);

% results = ["X","Y","Z","scalars"]
results = [msh.Points, A_total]
writeResults(results);

disp('Done!')