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

% plot mesh
figure(1)
tetramesh(msh);

[filename, filefolder] = uisavefile({'*.csv','Comma separated file';'*.txt','Raw Text File';'*.vtk','Binary Paraview file'}, 'Save output');
output_path = strcat(filefolder,filename);