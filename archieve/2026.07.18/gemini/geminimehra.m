clear; clc; close all;

%% 1. System Parameters & Grid Setup
R = 15;                             % Slew-rate limit [deg/s]
qco = 1;                            % Pilot command (kick) amplitude
fHzList   = 0.1:0.1:2;              % Command frequency [Hz]
wradsList = 2*pi*fHzList;
tauList   = [0, 0.03, 0.06, 0.09];  % Time delays [s]
KpList    = 0.1:0.1:5;              % Pilot gain grid

% Pre-allocate 3D array for storing Limit Cycle Amplitudes
% Dimensions: [Frequency, Pilot Gain, Time Delay]
ampGrid = zeros(length(wradsList), length(KpList), length(tauList));

%% 2. Prepare Simulink Model
mdl = "Mehra";
load_system(mdl);
set_param(mdl, "StopTime", "30");   % Sufficient time to reach steady-state
set_param(mdl, "SaveOutput", "on"); % Enable workspace output logging

disp('Running parametric sweep across w, Kp, and tau...');

%% 3. Nested Loop Simulation Sweep (Figure 7 Data)
for m = 1:length(tauList)
    tau = tauList(m);
    fprintf('Simulating for tau = %.2f s...\n', tau);
    
    for l = 1:length(KpList)
        Kp = KpList(l);
        
        for k = 1:length(wradsList)
            w = wradsList(k);
            
            % Push variables to base workspace for Simulink to read
            assignin('base', 'R', R);
            assignin('base', 'qco', qco);
            assignin('base', 'w', w);
            assignin('base', 'Kp', Kp);
            assignin('base', 'tau', tau);
            
            % Run simulation silently
            out = sim(mdl);
            
            % Robust signal extraction (handles different Simulink logging setups)
            if isprop(out, 'theta')
                theta_sig = out.theta.Data;
            elseif isprop(out, 'yout') && ~isempty(out.yout.find('theta'))
                theta_sig = out.yout.get('theta').Values.Data;
            elseif isprop(out, 'logsout') && ~isempty(out.logsout.find('theta'))
                theta_sig = out.logsout.get('theta').Values.Data;
            else
                error('Signal "theta" not found. Check Simulink logging instructions below.');
            end
            
            % Calculate steady-state amplitude using the last 40% of the time series
            idx_ss = round(length(theta_sig) * 0.6) : length(theta_sig);
            ampGrid(k, l, m) = (max(theta_sig(idx_ss)) - min(theta_sig(idx_ss))) / 2;
        end
    end
end

%% 4. Plot Figure 7: Jump Phenomena (Limit Cycle Amplitude vs. Kp)
% Extract the maximum amplitude across the frequency band for each (Kp, tau)
maxAmpOverFreq = squeeze(max(ampGrid, [], 1));

figure('Name', 'Figure 7: Jump Phenomena', 'Color', 'w', 'Position', [100, 100, 600, 500]);
hold on; grid on; box on;

% Match line styles from Mehra & Prasanth (1998)
lineStyles = {'--', '-.', ':', '-'};
legendLabels = cell(1, length(tauList));

for m = 1:length(tauList)
    plot(KpList, maxAmpOverFreq(:, m), 'LineWidth', 1.5, 'LineStyle', lineStyles{m});
    legendLabels{m} = sprintf('\\tau = %.2f', tauList(m));
end

xlabel('Pilot Gain K_p', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Limit Cycle Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
title('Figure 7: X-15 Limit cycle amplitude showing jump phenomena', 'FontSize', 11);
legend(legendLabels, 'Location', 'northwest', 'FontSize', 10);
xlim([0 5]);
ylim([0 max(maxAmpOverFreq(:))*1.1]);

%% 5. Plot Figure 8: Phase Portrait at Kp = 5.0
disp('Running dedicated simulation for Figure 8 (Phase Portrait)...');

% Set parameters for Kp = 5.0 limit cycle analysis
Kp = 5.0;
tau = 0;            % Baseline delay to generate clean spiral
w = 2*pi*0.5;       % Representative resonant forcing frequency
assignin('base', 'Kp', Kp);
assignin('base', 'tau', tau);
assignin('base', 'w', w);

% Extend stop time to capture full convergence from origin to limit cycle
set_param(mdl, "StopTime", "60");
out = sim(mdl);

% Extract theta and theta_dot
if isprop(out, 'yout') && ~isempty(out.yout.find('theta'))
    th     = out.yout.get('theta').Values.Data;
    th_dot = out.yout.get('theta_dot').Values.Data;
elseif isprop(out, 'logsout') && ~isempty(out.logsout.find('theta'))
    th     = out.logsout.get('theta').Values.Data;
    th_dot = out.logsout.get('theta_dot').Values.Data;
else
    th     = out.theta.Data;
    th_dot = out.theta_dot.Data;
end

figure('Name', 'Figure 8: Phase Portrait', 'Color', 'w', 'Position', [750, 100, 500, 500]);
plot(th, th_dot, 'k-', 'LineWidth', 1.0);
grid on; box on;
xlabel('\theta', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('d\theta/dt', 'FontSize', 12, 'FontWeight', 'bold');
title('Figure 8: X-15 Limit Cycle at K_p = 5.0', 'FontSize', 11);
axis equal;

% Add origin crosshairs
xline(0, 'k--', 'Alpha', 0.5);
yline(0, 'k--', 'Alpha', 0.5);