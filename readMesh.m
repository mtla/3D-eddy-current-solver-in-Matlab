function [ mesh ] = readMesh( source )
%READMESH Summary of this function goes here
%   Detailed explanation goes here
    if (ischar(source))
        mesh = 'string';
    elseif (ismatrix(source) && size(source,2) > 2)
        mesh = 'matrix';
    else
        mesh = 'error';
    end
end

