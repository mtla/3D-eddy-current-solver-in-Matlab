import matlab.unittest.TestSuite

% suite = TestSuite.fromClass(?tetrahedronTest);
% result = run(suite);
% 
% resultTable = table(result);
% sortrows(resultTable, 'Name')

suite = TestSuite.fromClass(?stiffnessMatrixTest);
result = run(suite);

resultTable = table(result);
sortrows(resultTable, 'Name')