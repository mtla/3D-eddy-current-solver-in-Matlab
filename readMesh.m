function [ mesh ] = readMesh( source )
%READMESH This function inputs vertices either in an .obj file, .txt file
%or an matrix/array. It then runs the mesh through a delaunay triangulation
%to split it into tetrahedrons. Finally, it returns the list of vertices
%and tetrahedrons in a struct.
%
% input:
%
%              nothing: opens a dialog so user can choose a file
%               string: tries to open that file
%         matrix/array: directly input a mesh as [n×3] matrix
%
% output:
%
%   delaunayTriangulation with properties:
%               Points: [n×3 double]
%     ConnectivityList: [m×4 double]
%          Constraints: [] (usually empty)
%
% More info of delaunayTriangulation here: https://se.mathworks.com/help/matlab/ref/delaunaytriangulation-class.html
%

    mesh = [];    
    if (exist('source','var')==0)
        [filename, filefolder] = uigetfile({'*.obj','Object Files';'*.txt','Raw Text File'}, 'Read obj-file');
        source = strcat(filefolder,filename);
    end
    if (ischar(source))
        [filefolder, filename, fileformat] = fileparts(source);
        fullname = (strcat(filefolder,'\',filename,fileformat));
        if (strcmp(fileformat,'.obj'))
            mesh = readObj(fullname);
        else % assume .txt file
            mesh = readTxt(fullname);
        end
    elseif (ismatrix(source) && size(source,2) > 2)
        mesh = source;
    end
    mesh = delaunayTriangulation(mesh);
end

