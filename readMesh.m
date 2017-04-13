function [ mesh ] = readMesh( source )
%READMESH This function inputs vertices either in an .obj file, .txt file
%or an matrix/array. It then runs the mesh through a delaunay triangulation
%to split it into tetrahedrons. Finally, it returns the list of vertices
%and tetrahedrons in a struct.
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
    mesh = delaunayTriangulation(mesh);
    
    
end

