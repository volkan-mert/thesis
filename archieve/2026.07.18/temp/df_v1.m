%% Nichols Plot Analysis: G(jw) and -1/N(A,w)
% Plots the frequency response of a linear plant against a family of 
% amplitude- and frequency-dependent describing function trajectories.
clear; clc; close all;

%% 1. Define Linear Plant G(jw)
% Example: Third-order plant G(s) = K / (s * (0.5*s + 1) * (0.1*s + 1))
K = 15;
w = logspace(-1, 2, 1000);          % Frequency vector (rad/s)
s = 1i * w;
G_jw = K ./ (s .* (0.5*s + 1) .* (0.1*s + 1));

% Convert Plant to Nichols Coordinates
mag_G_dB = 20 * log10(abs(G_jw));
phase_G_deg = rad2deg(angle(G_jw));

% Map phase continuously into the standard [-360, 0] Nichols window
phase_G_deg = mod(phase_G_deg, 360);
phase_G_deg(phase_G_deg > 0) = phase_G_deg(phase_G_deg > 0) - 360;

%% 2. Define Describing Function N(A, w)
% To demonstrate N(A,w), we use a Backlash nonlinearity (half-width 'b') 
% combined with an unmodeled time delay 'tau' (adds frequency-dependent phase lag).
b = 0.5;                            % Backlash dead-band half-width
tau = 0.04;                         % Time delay in seconds

% Backlash Describing Function N_b(A) valid for A > b
N_backlash = @(A) (1/pi) * (pi/2 + asin(1 - 2*b./A) + ...
    2*(1 - 2*b./A).*sqrt(b./A.*(1 - b./A))) ...
    - 1i * (4*b ./ (pi*A)) .* (1 - b./A);

% Total N(A, w) with frequency-dependent phase shift: N_b(A) * exp(-j*w*tau)
N_DF = @(A, w_val) N_backlash(A) .* exp(-1i * w_val * tau);

%% 3. Evaluate Family of -1/N(A, w) Trajectories
A_vec = linspace(b + 0.01, 5, 300); % Amplitude vector (> b)
w_samples = [1, 2, 3.5, 5, 7];      % Specific frequencies to evaluate N(A, w)

%% 4. Plotting on Nichols Chart
figure('Name', 'Nichols Plot: G(jw) vs -1/N(A,w)', 'Color', 'w', 'Position', [100, 100, 850, 600]);
hold on; grid on;

% Plot Linear Plant G(jw)
plot(phase_G_deg, mag_G_dB, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Plant G(j\omega)');

% Plot -1/N(A, w) family of curves for constant frequencies
colors = parula(length(w_samples) + 2);
for i = 1:length(w_samples)
    w_val = w_samples(i);

    % Compute negative inverse describing function
    neg_inv_N = -1 ./ N_DF(A_vec, w_val);

    % Convert to Nichols Coordinates (dB and Degrees)
    mag_N_dB = 20 * log10(abs(neg_inv_N));
    phase_N_deg = rad2deg(angle(neg_inv_N));

    % Ensure phase aligns near the -180 degree critical axis
    phase_N_deg = mod(phase_N_deg, 360);
    phase_N_deg(phase_N_deg > 0) = phase_N_deg(phase_N_deg > 0) - 360;

    % Plot trajectory for this specific frequency
    plot(phase_N_deg, mag_N_dB, '--', 'LineWidth', 1.8, 'Color', colors(i,:), ...
        'DisplayName', sprintf('-1/N(A, \\omega = %.1f rad/s)', w_val));

    % Mark the direction of increasing amplitude A with an arrow/marker at A = 1.5
    idx_marker = find(A_vec >= 1.5, 1);
    plot(phase_N_deg(idx_marker), mag_N_dB(idx_marker), 'o', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'HandleVisibility', 'off');
end

% Highlight the Critical Point (-180 deg, 0 dB)
plot(-180, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Critical Point (-180^\circ, 0 dB)');

%% 5. Chart Formatting & Grid Overlay
xlabel('Open-Loop Phase (degrees)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Open-Loop Magnitude (dB)', 'FontWeight', 'bold', 'FontSize', 11);
title('Nichols Plot: Harmonic Balance Analysis G(j\omega) = -1/N(A,\omega)', 'FontSize', 12);
axis([-300, -60, -30, 30]);         % Focus on the intersection region
legend('Location', 'northeast', 'FontSize', 10);

% Overlay standard Nichols grid if Control System Toolbox is available
if license('test', 'Control_Toolbox')
    ngrid; 
    set(findall(gcf,'String','0 dB'),'String',''); % Clean up overlapping grid labels if necessary
end
hold off;