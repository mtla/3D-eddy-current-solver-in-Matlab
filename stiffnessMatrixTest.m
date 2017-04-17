classdef stiffnessMatrixTest < matlab.unittest.TestCase
    % STIFFNESSMATRIXTEST
    % Tests whether stiffness matrixes are constructed correctly
    
    properties (TestParameter)
        meshes = struct(...
            'smallCube1', [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1], ...
            'smallCube2', [0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1], ...
            'randomCube1', randi(10, 8, 3),...
            'randomCube2', randi(10, 8, 3),...
            'randomCube3', randi(10, 8, 3),...
            'randomCube4', randi(10, 8, 3)...
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
        function testDimensions(testCase, meshes)
            DT = readMesh(meshes);
            actOutput = buildStiffnessMatrix(DT.ConnectivityList, DT.Points);
            actSolution = size(actOutput);
            expSolution = ones(1,2) * size(DT.Points, 1);
            testCase.verifyEqual(actSolution, expSolution, ... 
                'Return matrix should be an nxn matrix where n is the number of elements');
        end
    end
end

