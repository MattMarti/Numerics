function Test_All()
% Test_All.m
% 
% Recursively searches the Matlab folder in the Hybrid-Motor-Analysis repo
% and runs scripts that begin with the characters "test_" (Case
% insensitive).
%
% @authors: Matt Marti
% @date: 2019-08-31

clear, clc, close all

%% Recursively search for test scripts

% Determine repo location in the directory tree
m = inmem('-completenames');
full_testall_path = m{contains(m, 'Test_All.m')};
idum = regexp(full_testall_path, filesep);
path = full_testall_path(1:idum(end));

% Recursively run test scripts
recursive_test(path);


end

function recursive_test(str)
% Recursively add all directories in str to the path

% Find and run tests
file_list = dir(str);
for ii = 3:length(file_list)
    file = file_list(ii).name;
    if strncmpi('test_', file, 5) ...
            && strcmpi('.m', file(end-1:end)) ...
            && ~strcmp('Test_All.m', file)
        try
            run_test_script(file)
            fprintf('PASSED: %s\n', file(1:end-2))
        catch err
            fprintf('FAILED: %s\n\t%s\n', file(1:end-2), err.message)
        end
    end
end

% Perform action in folders
dirstruct = dir(str);
for i = 3:length(dirstruct)
    d = dirstruct(i);
    dirname = [str filesep d.name];
    if exist(dirname, 'dir')
        recursive_test( dirname );
    end
end

end

function run_test_script(file)

eval(file(1:end-2));

end
