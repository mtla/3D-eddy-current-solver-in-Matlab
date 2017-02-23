function vertices = readVerts( filename )
%READVERTS Reads verts from file with three columns
%   Input: text file with 3 columns separated by whitespace
%   Output: [nx3] matrix of the vertices
fid = fopen(filename);

tline = fgetl(fid);

CX = [];
CY = [];
CZ = [];

while ischar(tline)
%     % skip < and >
%     tline = substr(tline, 1, length(tline)-2)

    % extract numbers
    temp = textscan(tline,'%f%f%f');
%     celldisp(temp);
%     disp(temp(1));
    CX(end+1) = temp{1};
    CY(end+1) = temp{2};
    CZ(end+1) = temp{3};

    tline = fgetl(fid);
end

fclose(fid);
vertices = [CX' CY' CZ'];
end

