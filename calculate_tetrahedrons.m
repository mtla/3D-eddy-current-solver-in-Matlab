% function tetrahedrons = calculate_tetrahedrons()
% This function converts an array of 3D vertices into
% the appropriate formulation needed in the FEM calculation
%
% input: [m×3 double]
% output: [m×3 double]
%
% example input: 
%      1     1    -1
%      1    -1    -1
%      1     1     1
%      1    -1     1
%     -1     1    -1
%     -1    -1    -1
%     -1     1     1
%     -1    -1     1
%
% exmaple output:
% 
    

    verts = read_obj;
    v = verts.vertices;

    DT = delaunayTriangulation(v);
    tetrahedrons = DT;
    tetramesh(DT);
% end

