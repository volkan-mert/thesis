%% Rate Limiter Element (RLE) Nichols Chart via nicholsplot()
clear; clc; close all;

%% 1. System Parameters & Frequency Vector
R = 15;                       % Slew rate / Max actuator rate (deg/s)
A = 1;                        % Input amplitude (deg)
w = logspace(-1, 2, 1000);    % Frequency vector from 0.1 to 100 rad/s

% Calculate onset and critical saturation frequencies
w_onset = R / A;              % w_onset = 15 rad/s
w_crit  = 1.862 * w_onset;    % Fully developed saturation limit = 27.93 rad/s

%% 2. Describing Function Evaluation Across Regions
mag = zeros(1, length(w));
phi = zeros(1, length(w));

for k = 1:length(w)
    alpha = w(k) / w_onset;   % Frequency ratio (w / w_onset)

    if alpha < 1
        % Region I: No saturation (Eq. 15)
        mag(k) = 1;
        phi(k) = 0;

    elseif alpha < 1.862
        % Region II: Transition via cubic spline interpolation (Eq. 16)
        mag(k) = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
        phi(k) = 0.5280*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171; % rad

    else
        % Region III: Fully developed saturation (Eq. 17)
        varpi = 1 / alpha;    % varpi = w_onset / w
        mag(k) = (4 * varpi) / pi;
        phi(k) = -acos(pi * varpi / 2);                                 % rad
    end
end

%% 3. Package as Frequency Response Data (FRD) Model
% Combine linear magnitude and phase (in radians) into complex response array
N_jw = mag .* exp(1j * phi);

% Create the FRD LTI object (Requires Control System Toolbox)
sys_rle = frd(N_jw, w);

%% 4. Evaluate Boundary Markers for Visualization
% Saturation onset point (w = w_onset -> Region I boundary)
mag_onset_dB  = 20 * log10(1); 
phi_onset_deg = 0;

% Critical saturation point (w = 1.862 * w_onset -> Region III boundary)
alpha_c = 1.862;
mag_c   = 0.2908*alpha_c^3 - 1.4396*alpha_c^2 + 1.9232*alpha_c + 0.223;
phi_c   = 0.5280*alpha_c^3 - 2.6213*alpha_c^2 + 3.5056*alpha_c - 1.4171;
mag_crit_dB  = 20 * log10(mag_c);
phi_crit_deg = rad2deg(phi_c);

%% 5. Generate Nichols Chart using nicholsplot()
figure('Name', 'RLE Nichols Chart via nicholsplot', 'Color', 'w', 'Position', [100, 100, 850, 600]);

% Configure plot options natively for Control System Toolbox charts
opts = nicholsoptions('crossover');
opts.Title.String = 'Rate Limiter Element Describing Function N(j\omega, A_i)';
opts.Title.FontSize = 12;
opts.Title.FontWeight = 'bold';
opts.XLabel.FontSize = 11;
opts.YLabel.FontSize = 11;
opts.Grid = 'on';
opts.PhaseUnits = 'deg';
opts.MagUnits = 'dB';

% Plot the FRD LTI model
hPlot = nicholsplot(sys_rle, opts);

% Customize line appearance of the generated plot handle
findobj(gcf, 'Type', 'Line');
set(findall(gcf, 'Type', 'Line', 'DisplayName', ''), 'LineWidth', 2.5, 'Color', 'b');

% Overlay boundary markers directly onto the Nichols axes (Phase vs Gain)
hold on;
plot(phi_onset_deg, mag_onset_dB, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', ...
    'LineWidth', 1.5, 'DisplayName', '\omega_{onset} = 15 rad/s (Saturation Onset)');
plot(phi_crit_deg, mag_crit_dB, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
    'LineWidth', 1.5, 'DisplayName', '1.862\omega_{onset} = 27.93 rad/s (Full Saturation)');

% Add clear region text labels
text(phi_onset_deg - 3, mag_onset_dB + 1.2, 'Region I \rightarrow II', 'FontSize', 10, 'FontWeight', 'bold');
text(phi_crit_deg - 3, mag_crit_dB - 1.5, 'Region II \rightarrow III', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');

% Final axes formatting
axis([-100, 10, -25, 5]);
legend('RLE Describing Function N(j\omega)', '\omega_{onset} (15 rad/s)', '1.862\omega_{onset} (27.93 rad/s)', 'Location', 'southwest', 'FontSize', 10);