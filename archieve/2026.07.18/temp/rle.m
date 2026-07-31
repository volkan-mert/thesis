%% Rate Limiter Element (RLE) - Nichols Chart Analysis
clear; clc; close all;

%% 1. System Parameters & Frequency Vector
R = 15;                       % Slew rate / Max actuator rate (deg/s)
A = 1;                        % Input amplitude (deg)
w = logspace(-1, 2, 1000);    % Frequency vector from 0.1 to 100 rad/s

% Calculate onset and critical frequencies
w_onset = R / A;              % w_onset = 15 rad/s
w_crit  = 1.862 * w_onset;    % Fully developed saturation limit = 27.93 rad/s

%% 2. Describing Function Evaluation Across Regions
mag = zeros(size(w));
phi = zeros(size(w));

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
        omega_bar = 1 / alpha;    % omega_bar = w_onset / w
        mag(k) = (4 * omega_bar) / pi;
        phi(k) = -acos(pi * omega_bar / 2);                                 % rad
    end
end

% Convert to dB and degrees for Nichols presentation
mag_dB  = 20 * log10(mag);
phi_deg = rad2deg(phi);

%% 3. Evaluate Key Boundary Markers
% Onset point (w = w_onset)
mag_onset = 0; 
phi_onset = 0;

% Critical saturation point (w = 1.862 * w_onset)
alpha_crit = 1.862;
mag_crit_val = 0.2908*alpha_crit^3 - 1.4396*alpha_crit^2 + 1.9232*alpha_crit + 0.223;
phi_crit_val = 0.5280*alpha_crit^3 - 2.6213*alpha_crit^2 + 3.5056*alpha_crit - 1.4171;
mag_crit_dB  = 20 * log10(mag_crit_val);
phi_crit_deg = rad2deg(phi_crit_val);

%% 4. Plot Nichols Chart
figure('Name', 'RLE Nichols Chart', 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on
% Plot RLE describing function trajectory
plot(phi_deg, mag_dB, 'b', 'LineWidth', 2.5, 'DisplayName', 'RLE Describing Function N(j\omega, A_i)');

% Highlight Region Boundaries
plot(phi_onset, mag_onset, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', ...
    'DisplayName', '\omega_{onset} = 15 rad/s (Saturation Onset)');
plot(phi_crit_deg, mag_crit_dB, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
    'DisplayName', '1.862\omega_{onset} = 27.93 rad/s (Full Saturation)');

% Annotations & Formatting
text(phi_onset - 2, mag_onset + 1, 'Region I \rightarrow II (\omega_{onset})', 'FontSize', 10, 'FontWeight', 'bold');
text(phi_crit_deg - 2, mag_crit_dB - 1.5, 'Region II \rightarrow III (1.862\omega_{onset})', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');

title('Nichols Chart of Rate Limiter Element (RLE) Describing Function', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Open-Loop Phase (deg)', 'FontSize', 11);
ylabel('Open-Loop Gain (dB)', 'FontSize', 11);
xlim([-100, 10]);
ylim([-25, 5]);
legend('Location', 'southwest', 'FontSize', 10);

grid on

hold off