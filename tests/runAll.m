% Run all MATLAB tests in this folder.
testsFolder = fileparts(mfilename("fullpath"));
results = runtests(testsFolder);
disp(table(results));
assertSuccess(results);
