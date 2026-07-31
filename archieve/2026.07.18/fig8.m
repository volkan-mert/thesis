clear; clc; close all;

%% 1. Linear Aircraft Plant Definition
num_claw   = [5.21, -273.7855, -1425.456, -700.224];
den_claw   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw    = tf(num_claw, den_claw);

num_ldynac = [-10.524, -16.8384, -0.6209, 0];
den_ldynac = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac  = tf(num_ldynac, den_ldynac);

Gs_ac      = Gs_claw * Gs_ldynac;               % Linear OLTF q / q_c

%% 2. Rate Limiter Element (RLE) Describing Function Evaluation
R = 60;                       % Slew rate / Max actuator rate (deg/s)
A = 13.68;                        % Input amplitude (deg)
w = logspace(-1, 2, 1000);    % Frequency vector from 0.1 to 100 rad/s

% Calculate onset and critical saturation frequencies
w_onset = R / A;              % Saturation onset = 15 rad/s
w_crit  = 1.862 * w_onset;    % Fully developed saturation limit = 27.93 rad/s

mag = zeros(1, length(w));
phi = zeros(1, length(w));

for k = 1:length(w)
    alpha = w(k) / w_onset;   
    if alpha < 1
        % Region I: No saturation (Eq. 15)
        mag(k) = 1; phi(k) = 0;
    elseif alpha < 1.862
        % Region II: Transition via cubic spline interpolation (Eq. 16)
        mag(k) = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
        phi(k) = 0.5280*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171; % rad
    else
        % Region III: Fully developed saturation (Eq. 17)
        wbar = 1 / alpha;    
        mag(k) = (4 * wbar) / pi;
        phi(k) = -acos(pi * wbar / 2);                                  % rad
    end
end

%% 3. FRD Model Creation & Naming
N_jw = mag .* exp(1j * phi);
sys_rle = frd(N_jw, w);
indf_rle = -1/sys_rle; 


%% 4. Generate Aligned Nichols Chart via nicholsplot()
figure('Name', 'RLE Limit Cycle Analysis via nicholsplot', 'Color', 'w', 'Position', [100, 100, 900, 650]);

% Configure native Control System Toolbox plot options
opts = nicholsoptions('crossover');
opts.Title.String = 'Negative Inverse Describing Function Technique';
opts.Title.FontSize = 12;
opts.Title.FontWeight = 'bold';
opts.XLabel.FontSize = 11;
opts.YLabel.FontSize = 11;
opts.Grid = 'on';
opts.PhaseUnits = 'deg';
opts.MagUnits = 'dB';

% Force alignment of both models onto the 540 deg branch
opts.PhaseMatching = 'on';
opts.PhaseMatchingFreq = 15;      % Align around w_onset (15 rad/s)
opts.PhaseMatchingValue = 540;    % Target the 540 deg branch

% Plot native closed-loop contours and curves
nicholsplot(Gs_ac, indf_rle, opts);
legend('G(j\omega): Linear A/C TF','-1/N(j\omega, A_{i}): Negative Inverse DF of RLE')
hold on; 

%% 5. Exact PIO Limit Cycle Calculation (Without polyxpoly)
% Extract raw LTI frequency response of the aircraft
[mag_ac, phase_ac, w_ac] = nichols(Gs_ac, w);
mag_ac_dB = 20*log10(squeeze(mag_ac));
phase_ac_deg = squeeze(phase_ac);

% Negative inverse describing function (-1/N) shifted to 540° branch
mag_indf_dB = 20 * log10(1 ./ mag);
phase_indf_540 = (180 - rad2deg(phi)) + 360;  

% Use custom base-MATLAB intersection finder
[int_phase, int_mag, int_idx] = find_curve_intersect(phase_ac_deg, mag_ac_dB, phase_indf_540, mag_indf_dB);

%% 6. Overlay Frequency and Marker Directly on the Nichols Chart
if ~isempty(int_phase)
    % Interpolate frequency using exact fractional index along the plant curve
    idx_floor = floor(int_idx(1));
    t_frac    = int_idx(1) - idx_floor;
    w_pio     = w_ac(idx_floor) + t_frac * (w_ac(idx_floor+1) - w_ac(idx_floor));
    
    % Plot intersection marker
    plot(int_phase(1), int_mag(1), 'r*', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'LineWidth', 1.5, 'DisplayName', 'Limit Cycle Intersection');
    
    % Format and overlay semi-transparent text box over grid lines
    label_str = sprintf('  \\omega_{PIO} = %.2f rad/s (%.2f Hz)\n  Gain = %.2f dB\n  Phase = %.1f^\\circ', w_pio, w_pio/(2*pi), int_mag(1), int_phase(1));
                        
    text(int_phase(1), int_mag(1), label_str, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'k', ...
        'BackgroundColor', [1 1 1 0.85], ... % 85% opaque white box hides background grid
        'EdgeColor', 'k', ...                % Black border around text box
        'Margin', 5, ...
        'VerticalAlignment', 'bottom', ...   % Anchors box just above marker
        'HorizontalAlignment', 'left');
        
    % Print summary to Command Window
    fprintf('\n--- Predicted Limit Cycle (PIO) Point ---\n');
    fprintf('Oscillation Freq   : %.4f rad/s (%.4f Hz)\n', w_pio, w_pio/(2*pi));
    fprintf('Intersection Gain  : %.4f dB\n', int_mag(1));
    fprintf('Intersection Phase : %.4f deg\n\n', int_phase(1));
end

axis([360 720 -40 50]); % Zoom in around the 540 deg branch
legend('show', 'Location', 'best');

%% 7. Using Plot Function

figure;
plot(phase_ac_deg, mag_ac_dB, phase_indf_540, mag_indf_dB)
xlabel('Open-Loop Gain (deg)')
ylabel('Open-Loop Phase (dB)')
title('Negative Inverse Describing Function Technique')
grid on

%% LOCAL HELPER FUNCTION : INTERSECTION POINT OF THE CURVES %%
function [x_int, y_int, idx1] = find_curve_intersect(x1, y1, x2, y2)
    % Vectorized 2D line-segment intersection algorithm without Mapping Toolbox
    x_int = []; y_int = []; idx1 = [];
    
    % Ensure column vectors
    x1 = x1(:); y1 = y1(:); x2 = x2(:); y2 = y2(:);
    
    % Direction vectors of each segment (length N-1 and M-1)
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    % Starting points of each segment
    x1_start = x1(1:end-1); y1_start = y1(1:end-1);
    x2_start = x2(1:end-1); y2_start = y2(1:end-1);
    
    % Determinant matrix for all segment pairs ((N-1) x (M-1))
    DET = dx1 .* dy2' - dy1 .* dx2';
    
    % Coordinate differences between segment start points ((N-1) x (M-1))
    dX = x2_start' - x1_start;
    dY = y2_start' - y1_start;
    
    % Solve for parametric coordinates t (along curve 1) and u (along curve 2)
    T = (dX .* dy2' - dY .* dx2') ./ DET;
    U = (dX .* dy1  - dY .* dx1 ) ./ DET;
    
    % Segments intersect if 0 <= T <= 1 and 0 <= U <= 1 (ignore parallel lines)
    valid = (T >= 0) & (T <= 1) & (U >= 0) & (U <= 1) & (abs(DET) > 1e-10);
    
    if any(valid(:))
        [i, j] = find(valid);
        % Sort intersections by order of occurrence along curve 1
        [i, sort_order] = sort(i);
        j = j(sort_order);
        
        x_int = zeros(length(i), 1);
        y_int = zeros(length(i), 1);
        idx1  = zeros(length(i), 1);
        
        for k = 1:length(i)
            r = i(k); c = j(k);
            t_val = T(r, c);
            x_int(k) = x1(r) + t_val * dx1(r);
            y_int(k) = y1(r) + t_val * dy1(r);
            idx1(k)  = r + t_val; % Returns fractional index for frequency interpolation
        end
    end
end