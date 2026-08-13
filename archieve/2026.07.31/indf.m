clear; clc; close all

%% Define Linear System Components G(s)
s = tf('s');

% Control Law Transfer Function
num_claw = 5.21 * [1, -52.55, -273.6, -134.4];
den_claw = [1, 21.36, 545.6, 605.7, 0];
G_claw = tf(num_claw, den_claw);

% Longitudinal Dynamics Transfer Function
num_ac = -10.524 * [1, 1.6, 0.059, 0];
den_ac = [1, 2.35, -5.31, 0.184, -0.041];
G_ac = tf(num_ac, den_ac);

Kp = 1; % The Pilot Gain by taking the pilot model inside the loop

% Combined Open-Loop Transfer Function G(s)
G = Kp * G_claw * G_ac; 

%% Calculate -1/N Describing Function Loci
n = 1000;
w = logspace(-1, 2, n); 
R = 15; % Rate limit (R=15, u_rle=0.31)
u_rle = [0.2, 0.3, 0.31, 0.5, 1, 5]; % Different amplitudes of the input signal of R.L.E. (for R = 15 deg)
% u_rle = [1, 1.25, 1.5, 2]; % Different amplitudes of the input signal of R.L.E. (for R = 60 deg)

sys_inv_list = cell(1, length(u_rle));

for k = 1:length(u_rle)
    A_i = u_rle(k); % Amplitude of the input signal of the R.L.E
    
    % Describing Function N for the fully developed (Region III)
    varpi = R ./ (A_i * w);
    N = (4/pi) * varpi .* exp(-1i * acos((pi/2) * varpi));
    N((A_i * w) <= R) = 1; % Linear region
    
    % Inverse DF (-1/N)
    inv_N = -1 ./ N;
    
    % Reshape to 3D for FRD object (SISO format)
    inv_N_3D = reshape(inv_N, 1, 1, []);
    sys_inv_list{k} = frd(inv_N_3D, w);
end

%% Plot G(jw) and -1/N on the Nichols Chart
figure('Color', 'w');

% Plot G(s) along with all -1/N curves
nyquistplot(G, sys_inv_list{:}, w);

% Legend Formatting (Fixed arrayfun call)
legend_labels = [{'G(j\omega)'}, arrayfun(@(u) sprintf('-1/N (u_{rle} = %.4f)', u), u_rle, 'UniformOutput', false)];
legend(legend_labels, 'Location', 'northeast');

% Title with Slew Rate Information
title_slew_rate = sprintf('Slew Rate, R = %d deg/s', R);
title({'The Nichols Chart of a Rate Limiting Element ', 'Open-Loop Plant G(j\omega) vs. -1/N Loci', ['\color{red} (' title_slew_rate ')']});
