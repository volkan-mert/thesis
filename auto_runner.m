clear;
clc;

%% 1. Paths

msysBash = 'C:\msys64\usr\bin\bash.exe';

workDirWindows = 'C:\AUTO_work\ab_test';
workDirMSYS = '/c/AUTO_work/ab_test';

%% 2. Create working directory

if ~exist(workDirWindows, 'dir')
    mkdir(workDirWindows);
end

%% 3. Create AUTO CLUI script

autoScriptWindows = fullfile(workDirWindows, 'run_ab.auto');

fid = fopen(autoScriptWindows, 'w');

fprintf(fid, 'demo(''ab'')\n');
fprintf(fid, 'r = run(''ab'')\n');
fprintf(fid, 'save(r,''ab_matlab'')\n');

fclose(fid);

%% 4. Build command

bashCommand = [ ...
    'export PATH=/mingw64/bin:$PATH && ' ...
    'source /auto-07p/cmds/auto.env.sh && ' ...
    'cd ' workDirMSYS ' && ' ...
    'auto run_ab.auto'];

cmd = [ ...
    '"' msysBash '" -lc "' ...
    bashCommand ...
    '"'];

%% 5. Run AUTO

fprintf('\nRunning AUTO-07p from MATLAB...\n\n');

[status, output] = system(cmd);

fprintf('%s\n', output);

%% 6. Check result

if status == 0
    fprintf('\nAUTO-07p finished successfully.\n');
else
    fprintf('\nAUTO-07p returned an error.\n');
    fprintf('Status = %d\n', status);
end