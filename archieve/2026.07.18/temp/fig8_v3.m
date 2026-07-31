clear; clc; close all;

%% --------- Linear plant : control law × aircraft dynamics ----------------
num_claw   = [5.21, -273.7855, -1425.456, -700.224];
den_claw   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw    = tf(num_claw, den_claw);

num_ldynac = [-10.524, -16.8384, -0.6209, 0];
den_ldynac = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac  = tf(num_ldynac, den_ldynac);

Gs_ac      = Gs_claw * Gs_ldynac;          % linear OLTF q / q_c

%% Rate Limiter Element (RLE) Nichols Chart via nicholsplot()

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
        wbar = 1 / alpha;    % wbar = w_onset / w
        mag(k) = (4 * wbar) / pi;
        phi(k) = -acos(pi * wbar / 2);                                 % rad
    end
end

%% 3. Package as Frequency Response Data (FRD) Model
% Combine linear magnitude and phase (in radians) into complex response array
N_jw = mag .* exp(1j * phi);

% Create the FRD LTI object (Requires Control System Toolbox)
sys_rle = frd(N_jw, w);
indf_rle = -1/sys_rle; 

%% 4. Evaluate Boundary Markers for Visualization
mag_onset_dB  = 20 * log10(1); 
phi_onset_deg = 0;

alpha_c = 1.862;
mag_c   = 0.2908*alpha_c^3 - 1.4396*alpha_c^2 + 1.9232*alpha_c + 0.223;
phi_c   = 0.5280*alpha_c^3 - 2.6213*alpha_c^2 + 3.5056*alpha_c - 1.4171;
mag_crit_dB  = 20 * log10(mag_c);
phi_crit_deg = rad2deg(phi_c);

%% 5. Generate Aligned Nichols Chart using nicholsplot() & PhaseMatching
figure('Name', 'Aligned RLE Nichols Chart via nicholsplot', 'Color', 'w', 'Position', [100, 100, 850, 600]);

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

% --- MODIFICATION START: Force Phase Matching to 180° Branch ---
opts.PhaseMatching = 'on';
opts.PhaseMatchingFreq = 15;      % Align around w_onset (15 rad/s)
opts.PhaseMatchingValue = 180;    % Force both curves onto the 180 deg branch
% --- MODIFICATION END ---

% Plot both the linear plant and negative inverse describing function aligned
nicholsplot(Gs_ac, indf_rle, opts);
legend('Linear Plant G(j\omega)', '-1/N(j\omega, A_i)', 'Location', 'best');

%% 6. MODIFICATION: Explicit Shift to 540° Branch & Exact PIO Calculation
% To plot on the 540° branch (as seen in your original chart) and calculate 
% the exact limit cycle intersection point programmatically:

figure('Name', 'Exact PIO Intersection on 540° Branch', 'Color', 'w', 'Position', [150, 150, 850, 600]);
ngrid; hold on;

% Extract raw LTI frequency response of the aircraft
[mag_ac, phase_ac, w_ac] = nichols(Gs_ac, w);
mag_ac_dB = 20*log10(squeeze(mag_ac));
phase_ac_deg = squeeze(phase_ac);

% Calculate raw negative inverse describing function in dB and degrees
mag_indf_dB = 20 * log10(1 ./ mag);
phase_indf_deg = 180 - rad2deg(phi);  % Starts at 180° and moves positive

% Shift the describing function by +360° to align with the 540° loop
phase_indf_540 = phase_indf_deg + 360;

% Plot the curves on the 540° branch
plot(phase_ac_deg, mag_ac_dB, 'b', 'LineWidth', 1.5, 'DisplayName', 'Linear Plant G(j\omega)');
plot(phase_indf_540, mag_indf_dB, 'r', 'LineWidth', 1.5, 'DisplayName', '-1/N(j\omega, A_i) (540^\circ Branch)');

% Find exact intersection using polyxpoly
[int_phase, int_mag] = polyxpoly(phase_ac_deg, mag_ac_dB, phase_indf_540, mag_indf_dB);

if ~isempty(int_phase)
    % Plot the critical PIO intersection point
    plot(int_phase(1), int_mag(1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'LineWidth', 1.5, 'DisplayName', 'Limit Cycle Intersection');
    
    % Interpolate to find exact PIO oscillation frequency (w_pio)
    w_pio = interp1(phase_ac_deg, w_ac, int_phase(1));
    
    % Display limit cycle characteristics on chart
    txt = sprintf('  PIO Limit Cycle:\n  \\omega_{PIO} = %.2f rad/s\n  Gain = %.2f dB', w_pio, int_mag(1));
    text(int_phase(1), int_mag(1), txt, 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
    
    fprintf('\n--- Predicted Limit Cycle (PIO) Point ---\n');
    fprintf('Intersection Phase : %.2f deg\n', int_phase(1));
    fprintf('Intersection Gain  : %.2f dB\n', int_mag(1));
    fprintf('Oscillation Freq   : %.2f rad/s (%.2f Hz)\n\n', w_pio, w_pio/(2*pi));
end

title('Exact Rate Limiter Limit Cycle Intersection (540^\circ Branch)');
xlabel('Open-Loop Phase (deg)');
ylabel('Open-Loop Gain (dB)');
axis([360 720 -40 50]); % Zoom in around the 540° branch
legend('Location', 'southwest');
grid on;