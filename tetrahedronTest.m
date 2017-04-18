classdef tetrahedronTest < matlab.unittest.TestCase
    %TETRAHEDRONTEST 
    % Test class for testing the matrix that is returned from
    % tetrahedron2matrix function
    
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
        function testMappingToGlobal(testCase, tetrahedrons, vertices)
            actSolution = map2global(tetrahedrons, vertices);
            expSolution = ones(3);
            testCase.verifyEqual(actSolution, expSolution, ... 
                'The mapping of tetrahedron failed!');
        end
    end
end

