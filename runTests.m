import matlab.unittest.TestSuite

suite = TestSuite.fromFolder(pwd);
result = run(suite);

resultTable = table(result);
sortrows(resultTable, 'Name')