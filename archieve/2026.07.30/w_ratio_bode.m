% MATLAB Script to compute Amplitude (A) and Phase (phi) across regions
clear; clc; close all;

% 1. Define the range for the frequency ratio: r = w / w_onset
w_ratio = linspace(0.1, 3.0, 1000); 

% 2. Initialize output vectors
A   = zeros(size(w_ratio));
phi = zeros(size(w_ratio));

% 3. Region definitions using logical masks
reg1 = (w_ratio < 1);
reg2 = (w_ratio >= 1) & (w_ratio < 1.862);
reg3 = (w_ratio >= 1.862);

% --- Region I: No saturation ---
A(reg1)   = 1;
phi(reg1) = 0;

% --- Region II: Transition (Cubic spline interpolation) ---
% Here, alpha = w / w_onset
alpha = w_ratio(reg2);
A(reg2)   = 0.2908*alpha.^3 - 1.4396*alpha.^2 + 1.9232*alpha + 0.223;
phi(reg2) = 0.528*alpha.^3 - 2.6213*alpha.^2 + 3.5056*alpha - 1.4171;

% --- Region III: Fully developed saturation ---
% Here, varpi = w_onset / w
varpi = 1 ./ w_ratio(reg3);
A(reg3)   = (4 * varpi) / pi;
phi(reg3) = -acos((pi * varpi) / 2);

% 4. Visualization
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% Amplitude plot
subplot(2,1,1);
plot(w_ratio, A, 'b-', 'LineWidth', 2);
hold on;
xline(1.0, '--k', 'Region I / II');
xline(1.862, '--k', 'Region II / III');
grid on;
xlabel('\omega / \omega_{onset}', 'FontSize', 11);
ylabel('Amplitude (A)', 'FontSize', 11);
title('Amplitude A vs \omega / \omega_{onset}', 'FontSize', 12);

% Phase plot
subplot(2,1,2);
plot(w_ratio, phi, 'r-', 'LineWidth', 2);
hold on;
xline(1.0, '--k', 'Region I / II');
xline(1.862, '--k', 'Region II / III');
grid on;
xlabel('\omega / \omega_{onset}', 'FontSize', 11);
ylabel('Phase \phi (rad)', 'FontSize', 11);
title('Phase \phi vs \omega / \omega_{onset}', 'FontSize', 12);