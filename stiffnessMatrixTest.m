classdef stiffnessMatrixTest < matlab.unittest.TestCase
    %STIFFNESSMATRIXTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tetrahedrons
        node_list
    end
    properties (TestParameter)
    end
    
    methods(TestMethodSetup)
        function createBasicTetrahedrons(testCase)
            testCase.tetrahedrons = [8 6 5 2;7 8 5 2;7 4 8 2;7 1 4 2;7 5 1 2;7 3 4 1];
            testCase.node_list = [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1];
        end
    end
    
    methods (Test)
        function testDimensions(testCase)
            actOutput = buildStiffnessMatrix(testCase.tetrahedrons, testCase.node_list);
            actSolution = size(actOutput);
            expSolution = ones(1,2) * size(testCase.tetrahedrons, 1);
            testCase.verifyEqual(actSolution, expSolution, ... 
                'Return matrix should be an nxn matrix where n is the number of tetrahedrons');
        end
    end
end

