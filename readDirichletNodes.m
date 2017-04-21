function [ nodes ] = readDirichletNodes( source, DT )
% This function inputs vertices either in an .obj file, .txt file
% or an matrix/array. It then evaluates if the given vertices match the
% ones in the actual mesh and stores their order in a vector.
%
% input:
%
%              nothing: opens a dialog so user can choose a file
%               string: tries to open that file
%         matrix/array: directly input a mesh as [n×3] matrix
%
% output:
%
%               nodes: [n×1 double]
%
    
    if (exist('source','var')==0)
        [filename, filefolder] = uigetfile({'*.obj','Object Files';'*.txt','Raw Text File'}, 'Read obj-file');
        source = strcat(filefolder,filename);
    end
    if (ischar(source))
        [filefolder, filename, fileformat] = fileparts(source);
        fullname = (strcat(filefolder,'\',filename,fileformat));
        if (strcmp(fileformat,'.obj'))
            nodes_raw = readObj(fullname);
        else % assume .txt file
            nodes_raw = readTxt(fullname);
        end
        [node_exists, nodes] = ismember(nodes_raw, DT.Points,'rows');
        % check that all given dirichlet nodes exist in source mesh
        if any(~node_exists) 
            throw(MException('Input:NodeMissing','One or more of the given vertices do not exist in source mesh.'));
        end
    elseif (ismatrix(source) && size(source,2) > 0)
        % check that no index is larger than the number of vertices in the source mesh
        if max(source) > size(source,2) 
            throw(MException('Input:ArrayIndexOutOfBounds','Node index is larger than number of source nodes.'));
        end
        nodes = source;
    end
end

