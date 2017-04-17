classdef tetrahedronTest < matlab.unittest.TestCase
    %TETRAHEDRONTEST 
    % Test class for testing the matrix that is returned from
    % tetrahedron2matrix function
    % 
    % properties defines the variables used and TestMethodSetup assigns
    % values to those variables. The variables should be resetted after
    % each individual test.
    
    properties (TestParameter)
        vertices = struct('vertices', randi(10, 50,  3)); % mesh with 50 vertices
        tetrahedrons = struct(...
            'tetrahedron1', randperm(50,4),...
            'tetrahedron2', randperm(50, 4),...
            'tetrahedron3', randperm(50, 4),...
            'tetrahedron4', randperm(50, 4),...
            'tetrahedron5', randperm(50, 4),...
            'tetrahedron6', randperm(50, 4),...
            'tetrahedron7', randperm(50, 4),...
            'tetrahedron8', randperm(50, 4),...
            'tetrahedron9', randperm(50, 4)...
            );
    end
    
    methods (TestMethodSetup)
        function MethodSetup(testCase)
            orig = rng;
            testCase.addTeardown(@rng, orig)
            rng(42) % Answer to the Ultimate Question of Life, the Universe, and Everything.
        end
    end
    
%     properties
%         tetrahedron
%         node_list
%     end
%     properties (TestParameter)
%     end
%     
%     methods(TestMethodSetup)
%         function createBasicTetrahedron(testCase)
%             testCase.tetrahedron = [1 2 3 4];
%             testCase.node_list = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
%         end
%     end
    
    methods (Test)
        function testDimensions(testCase, vertices, tetrahedrons)
            actSolution = size(tetrahedron2matrix(tetrahedrons, vertices));
            expSolution = [4 4];
            testCase.verifyEqual(actSolution, expSolution, ... 
                'Return matrix should be 4x4');
        end
        function testSymmetry(testCase, vertices, tetrahedrons)
            actSolution = tetrahedron2matrix(tetrahedrons, vertices);
            expSolution = actSolution';
            testCase.verifyEqual(actSolution, expSolution, ... 
                'The returned matrix has to be symmetric');
        end
    end
end

