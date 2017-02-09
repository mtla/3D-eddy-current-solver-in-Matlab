function tetrahedrons = calculate_tetrahedrons( vertices )
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
    tetrahedrons = vertices;
end

