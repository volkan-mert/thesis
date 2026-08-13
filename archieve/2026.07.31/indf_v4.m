clear; clc; close all

%% 1. Define Linear System Components G(s)
s = tf('s');

% Control Law Transfer Function
num_claw = 5.21 * [1, -52.55, -273.6, -134.4];
den_claw = [1, 21.36, 545.6, 605.7, 0];
G_claw = tf(num_claw, den_claw);

% Longitudinal Dynamics Transfer Function
num_ac = -10.524 * [1, 1.6, 0.059, 0];
den_ac = [1, 2.35, -5.31, 0.184, -0.041];
G_ac = tf(num_ac, den_ac);

Kp = 1; % Pilot Gain
G = Kp * G_claw * G_ac; % Combined Open-Loop Transfer Function

R = 15; % Rate limit (deg/s)

%% 2. Solve for Exact Intersection Frequency (w_lc)
% Helper function to evaluate complex frequency response G(jw)
eval_G = @(w_val) squeeze(freqresp(G, w_val));

% Objective function: Phase_G(w) - Phase_invN(|G(jw)|) = 0
phase_residual = @(w_val) ...
    mod(rad2deg(angle(eval_G(w_val))), 360) - ...
    (180 + rad2deg(acos(min(1, (pi^2) / (8 * abs(eval_G(w_val)))))));

% Coarse grid search to locate initial root guess near 2-3 rad/s
w_grid = logspace(-1, 2, 2000);
res_grid = arrayfun(phase_residual, w_grid);
idx_sign_change = find(diff(sign(res_grid)) ~= 0 & abs(res_grid(1:end-1)) < 30, 1, 'last');

if isempty(idx_sign_change)
    error('No intersection point found between G(jw) and -1/N locus.');
end

% Refine exact frequency using fzero
w_lc = fzero(phase_residual, w_grid(idx_sign_change));
f_lc = w_lc / (2 * pi);

% Calculate Exact Coordinates at Intersection
G_lc = eval_G(w_lc);
gain_lc_db = 20 * log10(abs(G_lc));
phase_lc_deg = mod(rad2deg(angle(G_lc)), 360);

% Calculate Limit Cycle Saturation Parameter and Input Amplitude u_rle
x_lc = (4 / pi) * abs(G_lc);
u_rle_lc = (R * x_lc) / w_lc;

% Display Results in Command Window
fprintf('\n================== EXACT INTERSECTION RESULTS ==================\n');
fprintf('Limit Cycle Frequency (w_lc) : %.4f rad/s (%.4f Hz)\n', w_lc, f_lc);
fprintf('Open-Loop Gain at Crossing  : %.2f dB\n', gain_lc_db);
fprintf('Open-Loop Phase at Crossing : %.2f deg\n', phase_lc_deg);
fprintf('Required Input Amplitude    : u_rle = %.4f deg\n', u_rle_lc);
fprintf('================================================================\n\n');

%% 3. Plot Nichols Chart and Overlay Intersection Point
n = 1000;
w = logspace(-1, 2, n); 
u_rle_vec = [0.2, 0.3, 0.31, 0.5, 1, 5, u_rle_lc, 50]; 

sys_inv_list = cell(1, length(u_rle_vec));
for k = 1:length(u_rle_vec)
    A_i = u_rle_vec(k);
    varpi = R ./ (A_i * w);
    arg_clamp = min(1, (pi/2) * varpi);
    N = (4/pi) * varpi .* exp(-1i * acos(arg_clamp));
    N((A_i * w) <= R) = 1;
    inv_N_3D = reshape(-1 ./ N, 1, 1, []);
    sys_inv_list{k} = frd(inv_N_3D, w);
end

figure('Color', 'w');
[~, phase_start] = bode(G, w(1)); 

opt = nicholsoptions;
opt.Grid = 'on';
opt.PhaseMatching = 'on';
opt.PhaseMatchingFreq = w(1);
opt.PhaseMatchingValue = phase_start - 270;

h = nicholsplot(G, sys_inv_list{:}, w);
setoptions(h, opt);
hold on;

% Highlight Intersection Point with Red Circle and Marker
plot(phase_lc_deg, gain_lc_db, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(phase_lc_deg, gain_lc_db, 'r+', 'MarkerSize', 10, 'LineWidth', 2);

% Annotate Plot
text_str = sprintf('  \\leftarrow Intersection Point\n  \\omega_{lc} = %.3f rad/s (%.3f Hz)\n  Gain = %.2f dB, Phase = %.1f^o\n  u_{rle, lc} = %.4f', ...
    w_lc, f_lc, gain_lc_db, phase_lc_deg, u_rle_lc);
text(phse_lc_deg, gain_lc_db, text_str, 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 10);

legend_labels = [{'G(j\omega)'}, arrayfun(@(u) sprintf('-1/N (u_{rle} = %.4f)', u), u_rle_vec, 'UniformOutput', false)];
legend(legend_labels, 'Location', 'northeast');

title_slew_rate = sprintf('Slew Rate, R = %d deg/s', R);
title({'The Nichols Chart of a Rate Limiting Element', ...
       'Open-Loop Plant G(j\omega) vs. -1/N Loci', ...
       ['\color{red} (' title_slew_rate ')']});