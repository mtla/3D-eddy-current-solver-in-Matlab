[filename, filefolder] = uigetfile({'*.obj','Object Files';'*.txt','Raw Text File'}, 'Read obj-file');
mesh_path = strcat(filefolder,filename);

[filename, filefolder] = uigetfile({'*.obj','Object Files';'*.txt','Raw Text File'}, 'Read obj-file');
dirichlet_path = strcat(filefolder,filename);

msh = readMesh(mesh_path);
dirichletNodes = readDirichletNodes(msh, dirichlet_path);

np = size(msh.Points,1);
if (exist('passes','var')==0)
    passes = 2;
end
if (exist('permittivity','var')==0)
    permittivity = 10;
end
if (exist('permeability','var')==0)
    permeability = 10;
end

prompt = {'Subdivide mesh # times: ','Permittivity in element:','Permeability in element:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {num2str(passes), num2str(permittivity), num2str(permeability)};
constants = inputdlg(prompt,dlg_title,num_lines,defaultans);
[passes, permittivity, permeability] = constants{:};

passes = str2num(passes);
msh.refine(passes);

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

msh.setPointValues(A_total);
msh.setEdgeValues(Aedges);

figure(2)
scatter3(msh, A_total);
plot3(msh, Aedges);

% results = ["X","Y","Z","scalars"]
writeResults(msh);

disp('Done!')