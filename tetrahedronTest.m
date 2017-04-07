classdef tetrahedronTest < matlab.unittest.TestCase
    %TETRAHEDRONTEST 
    % Test class for testing the matrix that is returned from
    % tetrahedron2matrix function
    % 
    % properties defines the variables used and TestMethodSetup assigns
    % values to those variables. The variables should be resetted after
    % each individual test.
    
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
        function testSymmetry(testCase)
            actSolution = tetrahedron2matrix(testCase.tetrahedron, testCase.node_list);
            expSolution = actSolution';
            testCase.verifyEqual(actSolution, expSolution, ... 
                'The returned matrix has to be symmetric');
        end
    end
end

