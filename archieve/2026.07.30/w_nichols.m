clear; clc; close all
%%

% 1. Frequency ratio vector definition
n = 1000;
w = logspace(-1, 2, n);  % w = w / w_onset (0.1 to 100)

% 2. Preallocate vectors
A   = zeros(size(w));
phi = zeros(size(w));   % in radians

% 3. Region Masks
reg1 = (w < 1.0);
reg2 = (w >= 1.0) & (w < 1.862);
reg3 = (w >= 1.862);

% --- Region I: No saturation ---
A(reg1)   = 1.0;
phi(reg1) = 0.0;

% --- Region II: Transition (Cubic spline) ---
alpha = w(reg2);
A(reg2)   = 0.2908*alpha.^3 - 1.4396*alpha.^2 + 1.9232*alpha + 0.223;
phi(reg2) = 0.5280*alpha.^3 - 2.6213*alpha.^2 + 3.5056*alpha - 1.4171;

% --- Region III: Fully developed saturation ---
varpi = 1.0 ./ w(reg3);
A(reg3)   = (4.0 * varpi) / pi;
phi(reg3) = -acos((pi * varpi) / 2.0);

% 4. Convert to dB and Degrees for Nichols Plot
mag_db  = 20 * log10(A);     % Gain in dB
phi_deg = phi * (180 / pi);  % Phase in degrees

% =========================================================================
% METHOD A: Direct Nichols Plotting with M- and N-circles grid (ngrid)
% =========================================================================
figure('Name', 'Nichols Chart', 'Color', 'w', 'Position', [150, 150, 800, 600]);

% Plot Nichols Curve
plot(phi_deg, mag_db, 'b-', 'LineWidth', 2.2, 'DisplayName', 'Describing Function N(w)');
hold on;

% Overlay Nichols Grid (Constant closed-loop magnitude and phase curves)
ngrid;

% Highlight region transition boundaries
idx_r1_r2 = find(w >= 1.0, 1, 'first');
idx_r2_r3 = find(w >= 1.862, 1, 'first');

plot(phi_deg(idx_r1_r2), mag_db(idx_r1_r2), 'ro', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Region I/II Boundary (\omega/\omega_{onset} = 1.0)');

plot(phi_deg(idx_r2_r3), mag_db(idx_r2_r3), 'go', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Region II/III Boundary (\omega/\omega_{onset} = 1.862)');

% Axis Limits and Labels
xlim([-100, 10]);
ylim([-40, 5]);
xlabel('Phase (degrees)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Magnitude (dB)', 'FontSize', 11, 'FontWeight', 'bold');
title('Nichols Chart for Actuator Saturation Describing Function', 'FontSize', 12);
legend('Location', 'southwest');
hold off;

% =========================================================================
% METHOD B: Using MATLAB Control System Toolbox (FRD Object)
% =========================================================================
% Complex describing function N(w) = A * exp(j*phi)
N_complex = A .* exp(1i * phi);
sys_frd   = frd(N_complex, w);

figure('Name', 'Nichols Chart', 'Color', 'w');

figure;
nichols(-1/sys_frd,'k--');
ngrid;
title('Nichols Chart of the Negative Inverse Describing Function (by using FRD)');