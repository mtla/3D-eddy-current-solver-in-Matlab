function [ mesh ] = readMesh( source )
%READMESH Summary of this function goes here
%   Detailed explanation goes here
    mesh = [];
    if (ischar(source))
        fileformat = split(source, '.');
        fileformat = fileformat(end);
        if (strcmp(fileformat,'obj'))
            mesh = readObj(source);
        elseif (strcmp(fileformat, 'txt'))
            mesh = readTxt(source);
        end
    elseif (ismatrix(source) && size(source,2) > 2)
        mesh = source;
    end
end

