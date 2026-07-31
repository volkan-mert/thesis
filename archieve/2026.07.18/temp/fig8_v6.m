clear; clc; close all;

%% 1. Linear Plant Definition
num_claw   = [5.21, -273.7855, -1425.456, -700.224];
den_claw   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw    = tf(num_claw, den_claw);

num_ldynac = [-10.524, -16.8384, -0.6209, 0];
den_ldynac = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac  = tf(num_ldynac, den_ldynac);

Gs_ac      = Gs_claw * Gs_ldynac;          % Linear OLTF q / q_c
Gs_ac.Name = 'Linear Plant G(j\omega)';    % Name for legend

%% 2. Rate Limiter Describing Function Evaluation
R = 15;                       % Slew rate / Max actuator rate (deg/s)
A = 1;                        % Input amplitude (deg)
w = logspace(-1, 2, 1000);    % Frequency vector from 0.1 to 100 rad/s
w_onset = R / A;              % Saturation onset = 15 rad/s

mag = zeros(1, length(w));
phi = zeros(1, length(w));

for k = 1:length(w)
    alpha = w(k) / w_onset;   
    if alpha < 1
        mag(k) = 1; phi(k) = 0;
    elseif alpha < 1.862
        mag(k) = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
        phi(k) = 0.5280*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171; 
    else
        wbar = 1 / alpha;    
        mag(k) = (4 * wbar) / pi;
        phi(k) = -acos(pi * wbar / 2);                                 
    end
end

%% 3. Package as Frequency Response Data (FRD) Model
N_jw = mag .* exp(1j * phi);
sys_rle = frd(N_jw, w);
indf_rle = -1/sys_rle; 
indf_rle.Name = '-1/N(j\omega, A_i)';       % Name for legend

%% 4. Generate Chart Using nicholsplot()
figure('Name', 'Nichols Chart with Intersection Frequency', 'Color', 'w', 'Position', [150, 100, 900, 650]);

% Configure native Control System Toolbox plot options
opts = nicholsoptions('crossover');
opts.Title.String = 'Rate Limiter Describing Function Intersection on Nichols Chart';
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

% Plot using nicholsplot
nicholsplot(Gs_ac, indf_rle, opts);
hold on; % Hold axes to overlay intersection marker and text box

%% 5. Extract Raw Data & Calculate Intersection Numerically
[mag_ac, phase_ac, w_ac] = nichols(Gs_ac, w);
mag_ac_dB = 20*log10(squeeze(mag_ac));
phase_ac_deg = squeeze(phase_ac);

% Negative inverse describing function (-1/N) shifted to 540° branch
mag_indf_dB = 20 * log10(1 ./ mag);
phase_indf_540 = (180 - rad2deg(phi)) + 360;  

% Find intersection without polyxpoly
[int_phase, int_mag, int_idx] = find_curve_intersect(phase_ac_deg, mag_ac_dB, phase_indf_540, mag_indf_dB);

%% 6. Overlay Frequency and Label on the nicholsplot
if ~isempty(int_phase)
    % --- FREQUENCY INTERPOLATION ---
    idx_floor   = floor(int_idx(1));
    t_frac      = int_idx(1) - idx_floor;
    w_intersect = w_ac(idx_floor) + t_frac * (w_ac(idx_floor+1) - w_ac(idx_floor));
    
    % --- PLOT MARKER ON NICHOLS CHART ---
    plot(int_phase(1), int_mag(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'y', ...
        'LineWidth', 1.5, 'DisplayName', 'Limit Cycle Intersection');
    
    % --- FORMAT AND WRITE TEXT BOX OVER GRID ---
    label_str = sprintf('  \\omega_{int} = %.2f rad/s (%.2f Hz)\n  Gain = %.2f dB\n  Phase = %.1f^\\circ', ...
                        w_intersect, w_intersect/(2*pi), int_mag(1), int_phase(1));
                        
    text(int_phase(1), int_mag(1), label_str, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'k', ...
        'BackgroundColor', [1 1 1 0.85], ... % 85% opaque white box hides background grid
        'EdgeColor', 'k', ...                % Black border
        'Margin', 4, ...
        'VerticalAlignment', 'bottom', ...   % Anchors text box just above the marker
        'HorizontalAlignment', 'left');
        
    % Print summary to Command Window
    fprintf('--- Limit Cycle Intersection (PIO) ---\n');
    fprintf('Frequency : %.4f rad/s (%.4f Hz)\n', w_intersect, w_intersect/(2*pi));
    fprintf('Gain      : %.4f dB\n', int_mag(1));
    fprintf('Phase     : %.4f deg\n\n', int_phase(1));
end

axis([360 720 -40 50]); % Focus view on the 540 deg branch
legend('show', 'Location', 'southwest');

%% LOCAL HELPER FUNCTION : Intersection of the curves %%
function [x_int, y_int, idx1] = find_curve_intersect(x1, y1, x2, y2)
    % Vectorized 2D line-segment intersection algorithm without Mapping Toolbox
    x_int = []; y_int = []; idx1 = [];
    x1 = x1(:); y1 = y1(:); x2 = x2(:); y2 = y2(:);
    
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    x1_start = x1(1:end-1); y1_start = y1(1:end-1);
    x2_start = x2(1:end-1); y2_start = y2(1:end-1);
    
    DET = dx1 .* dy2' - dy1 .* dx2';
    dX  = x2_start' - x1_start;
    dY  = y2_start' - y1_start;
    
    T = (dX .* dy2' - dY .* dx2') ./ DET;
    U = (dX .* dy1  - dY .* dx1 ) ./ DET;
    
    valid = (T >= 0) & (T <= 1) & (U >= 0) & (U <= 1) & (abs(DET) > 1e-10);
    
    if any(valid(:))
        [i, j] = find(valid);
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
            idx1(k)  = r + t_val; % Returns fractional segment index
        end
    end
end