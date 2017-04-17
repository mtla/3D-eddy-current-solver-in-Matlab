classdef geometryTest < matlab.unittest.TestCase
    % GEOMETRYTEST
    % Test class for testing if program can handle different kinds of input
    % meshes
    
%     properties
%         mesh
%         mesh2
%     end
    
    properties (ClassSetupParameter)
        generator = {'twister','combRecursive','multFibonacci'};
    end
    
    properties (MethodSetupParameter)
        seed = {0, 123, 4294967295};
    end
    
    properties (TestParameter)
        meshes = struct('smallCube', [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1], ...
            'cubeWithDuplicates', [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;1 1 1;-1 -1 1], ...
            'randomCube', randi(10, 3, 8));
        files = struct('obj', strcat(pwd,'\meshes\example_mesh_3D.obj'), ...
            'txt', strcat(pwd,'\meshes\cube_nodes.txt'), ...
            'gibberish', 'huehue', ...
            'gibberish2', 'blaoeuchaoeurcoha', ...
            'gibberishDots', 'oeuouou.hue', ...
            'gibberishMultipleDots', 'orocehuo.hue.hue.hue.hue', ...
            'gibberishSlashes', 'aoeuoeuo\oeuoeu/ooeuo', ...
            'gibberishDotsSlashes', 'oecuhaoeu/oeuoceuh/oeuch.hue');
    end

    methods (TestClassSetup)
        function ClassSetup(testCase, generator)
            orig = rng;
            testCase.addTeardown(@rng, orig)
            rng(0, generator)
        end
    end
    
    methods (TestMethodSetup)
        function MethodSetup(testCase, seed)
            orig = rng;
            testCase.addTeardown(@rng, orig)
            rng(seed)
        end
    end
    
%     methods(TestMethodSetup)
%         function createBasicTetrahedron(testCase)
%             testCase.mesh = [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;1 1 1;-1 -1 1];
%             testCase.mesh2= [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1];
%         end
%     end
    
    methods (Test)
        % this test gives duplicate rows on purpose to the readMesh
        % function
        function testDuplicates(testCase, meshes)
            actSolution = readMesh(meshes);
            actSolution = actSolution.Points;
            expSolution = unique(actSolution, 'rows', 'stable');
            testCase.verifyEqual(actSolution, expSolution, ... 
                "Mesh can't have duplicate vertices");
        end
        function testTxt(testCase, files)
            actSolution = readMesh(files);
            expSolution = [0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1];
            testCase.verifyEqual(actSolution.Points, expSolution, ... 
                ".txt input failed!");
        end
%         function testObj(testCase)
%             actSolution = readMesh(strcat(pwd,'\meshes\example_mesh_3D.obj'));
%             expSolution = testCase.mesh2;
%             testCase.verifyEqual(actSolution.Points, expSolution, ... 
%                 ".obj input failed!");
%         end
        function testMatrix(testCase)
            actSolution = readMesh([1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1]);
            expSolution = [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1];
            testCase.verifyEqual(actSolution.Points, expSolution, ... 
                "Input of matrix/array failed!");
        end
    end
end