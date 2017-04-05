classdef tetrahedronTest < matlab.unittest.TestCase
    %TETRAHEDRONTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tetrahedron
        node_list
    end
    properties (TestParameter)
    end
    
    methods(TestMethodSetup)
        function createBasicTetrahedron(testCase)
            testCase.tetrahedron = [1 2 3 4];
            testCase.node_list = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
        end
    end
    
    methods (Test)
        function testDeterminant(testCase)
            actSolution = det(tetrahedron2matrix(testCase.tetrahedron, testCase.node_list));
            expSolution = 0;
            testCase.verifyEqual(actSolution, expSolution, ... 
                'Determinant of tetrahedron matrix should be 0');
        end
    end
    
end

