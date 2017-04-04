classdef tetrahedronTest < matlab.unittest.TestCase
    %TETRAHEDRONTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tetrahedron = [1 2 3 4];
        node_list = [0 0 0, 0 0 1, 0 1 0, 1 0 0];
    end
    
    methods (Test)
        function testDeterminant(testCase)
            actSolution = det(tetrahedron2matrix(tetrahedron, node_list));
            exSolution = 0;
            testCase.verifyEqual(actSolution, expSolution);
        end
    end
    
end

