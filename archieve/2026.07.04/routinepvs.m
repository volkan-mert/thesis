clear; clc ;close all

% 1. Define the model name and parameters
modelName = 'pvsa';
n = 10000;

% Create the logarithmically spaced frequency vector (0.1 to 100 rad/s)
w_values = logspace(-1, 2, n); 

% 2. Load the model into memory without opening the UI
load_system(modelName);

% 3. Preallocate the SimulationInput array for memory efficiency
% (Preallocating is critical when n = 10000 to prevent MATLAB from slowing down)
simIn(1:n) = Simulink.SimulationInput(modelName);

% 4. Populate the array with your varying 'w' parameters
for i = 1:n
    % Overwrite the variable 'w' in the workspace for this specific run
    simIn(i) = simIn(i).setVariable('w', w_values(i));
    
    % CRITICAL FOR PERFORMANCE: Enable Fast Restart
    % This skips the compile phase for runs 2 through 10,000
    simIn(i) = simIn(i).setModelParameter('FastRestart', 'on');
end

% 5. Execute the simulations in parallel
% Note: This requires the Parallel Computing Toolbox
fprintf('Starting parallel execution of %d simulations...\n', n);
out = parsim(simIn, 'ShowProgress', 'on');

% 6. (Optional) Extracting the data
% Assuming you have signal logging enabled for 'theta' or 'q'
% preallocate a cell array to extract results quickly:

theta_responses = cell(n, 1);
for i = 1:n
    % Check if the run was successful before extracting
    if empty(out(i).ErrorMessage)
        theta_responses{i} = out(i).logsout.get('theta').Values;
    else
        warning('Run %d failed: %s', i, out(i).ErrorMessage);
    end
end