classdef geometryTest < matlab.unittest.TestCase
    % GEOMETRYTEST
    % Test class for testing if program can handle different kinds of input
    % meshes
    
    properties
        mesh
        mesh2
    end
%     properties (ClassSetupParameter)
%         generator = {'twister','combRecursive','multFibonacci'};
%     end
%     
%     properties (MethodSetupParameter)
%         seed = {0, 123, 4294967295};
%     end
%     
%     properties (TestParameter)
%         dim1 = struct('small', 1,'medium', 2, 'large', 3);
%         dim2 = struct('small', 2,'medium', 3, 'large', 4);
%         dim3 = struct('small', 3,'medium', 4, 'large', 5);
%         type = {'single','double'};
%     end
    
    methods(TestMethodSetup)
        function createBasicTetrahedron(testCase)
            testCase.mesh = [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;1 1 1;-1 -1 1];
            testCase.mesh2= [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1];
        end
    end
    
    methods (Test)
        % this test gives duplicate rows on purpose to the readMesh
        % function
        function testDuplicates(testCase)
            actSolution = readMesh(testCase.mesh);
            expSolution = unique(actSolution, 'rows');
            testCase.verifyEqual(actSolution, expSolution, ... 
                "Mesh can't have duplicate vertices");
        end
        function testTxt(testCase)
            actSolution = readMesh(strcat(pwd,'\meshes\cube_nodes.txt'));
            expSolution = [0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1];
            testCase.verifyEqual(actSolution, expSolution, ... 
                ".txt input failed!");
        end
        function testObj(testCase)
            actSolution = readMesh(strcat(pwd,'\meshes\example_mesh_3D.obj'));
            expSolution = testCase.mesh2;
            testCase.verifyEqual(actSolution, expSolution, ... 
                ".obj input failed!");
        end
        function testMatrix(testCase)
            actSolution = readMesh([1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1]);
            expSolution = [1 1 -1;1 -1 -1;1 1 1;1 -1 1;-1 1 -1;-1 -1 -1;-1 1 1;-1 -1 1];
            testCase.verifyEqual(actSolution, expSolution, ... 
                "Input of matrix/array failed!");
        end
    end
end