classdef geometryTest < matlab.unittest.TestCase
    % GEOMETRYTEST
    % Test class for testing if program can handle different kinds of input
    % meshes
    
    properties
        mesh
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
        end
    end
    
    methods (Test)
        function testDuplicates(testCase)
            actSolution = readMesh(testCase.mesh);
            expSolution = unique(actSolution, 'rows');
            testCase.verifyEqual(actSolution, expSolution, ... 
                "Mesh can't have duplicate vertices");
        end
    end
end