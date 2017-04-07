classdef stiffnessMatrixTest < matlab.unittest.TestCase
    %STIFFNESSMATRIXTEST Summary of this class goes here
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
        function testDimensions(testCase)
            actSolution = size(tetrahedron2matrix(testCase.tetrahedron, testCase.node_list));
            expSolution = [4 4];
            testCase.verifyEqual(actSolution, expSolution, ... 
                'Return matrix should be 4x4');
        end
    end
end

