import matlab.unittest.TestSuite

suite = TestSuite.fromClass(?tetrahedronTest);
result = run(suite);

resultTable = table(result);
sortrows(resultTable, 'Name')