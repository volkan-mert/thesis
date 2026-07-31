%% APPLIED NONLINEAR CONTROL: DESCRIBING FUNCTION ANALYSIS
% Frequency-Dependent Describing Function (FDDF) Limit Cycle Detection
% Plant Model: Gs_ac = Gs_ldynac

clear; clc; close all;

%% 1. Linear Aircraft Plant Definition (Gs_ldynac only)
num_ldynac = [-10.524, -16.8384, -0.6209, 0];
den_ldynac = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac  = tf(num_ldynac, den_ldynac);

Gs_ac      = Gs_ldynac;               % Linear OLTF set to Gs_ldynac

%% 2. Rate Limiter Element (RLE) Describing Function Evaluation
R = 60;                       % Slew rate / Max actuator rate (deg/s)
A = 13.68;                    % Input amplitude (deg)
w = logspace(-1, 2, 2000);    % Frequency vector (rad/s)

% Calculate onset and critical saturation frequencies
w_onset = R / A;              
w_crit  = 1.862 * w_onset;    

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

%% 3. FRD Model Creation & Plotting
N_jw = mag .* exp(1j * phi);
sys_rle = frd(N_jw, w);
indf_rle = -1/sys_rle; 

figure;
nyquistplot(Gs_ac, indf_rle);
title('FREQUENCY-DEPENDENT DESCRIBING FUNCTION (Gs\_ldynac)', 'FontWeight', 'bold');
grid on;

%% 4. Mark the critical point (-1, 0)
hold on;
plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
text(-1.05, 0.2, ' (-1, 0i)', 'Color', 'r', 'FontWeight', 'bold');
ylim([-5, 5]);
xlim([-5, 5]);

%% 5. Geometric Intersection Search (Enforcing w_ac == w_rle)

% Extract complex frequency response arrays
Gs_ac_vec    = squeeze(freqresp(Gs_ac, w));
indf_rle_vec = squeeze(freqresp(indf_rle, w));

% Compute distance ONLY at matching frequencies across vector elements
dist_vec = abs(Gs_ac_vec - indf_rle_vec);

% Find frequency index with global minimum difference
[min_dist, best_k] = min(dist_vec);

% Distance tolerance threshold
tol = 0.15; 
HAS_INTERSECTION = (min_dist < tol);

if HAS_INTERSECTION
    % Coordinates and frequency at limit cycle condition
    x_int = real(indf_rle_vec(best_k));
    y_int = imag(indf_rle_vec(best_k));
    w_osc = w(best_k);
    A_rle = A;
    A_ac  = abs(Gs_ac_vec(best_k));

    % Plot red star at intersection
    plot(x_int, y_int, 'r*', 'MarkerSize', 12, 'LineWidth', 2);

    % Text callout for w_osc and A_rle
    text_str = sprintf('  \\omega_{osc} = %.2f rad/s\n  A_{rle} = %.2f deg', w_osc, A_rle);
    text(x_int + 0.2, y_int - 0.3, text_str, 'FontSize', 9, 'FontWeight', 'bold');

    % Plot Legend
    legend('Aircraft Transfer Function: Gs\_ldynac(j\omega)', ...
           'Inverse Negative Describing Function: (-1 / N(A, j\omega))', ...
           'Critical Point (-1, 0i)', ...
           'Limit Cycle Point', ...
           'Location', 'best');
else
    % Text box displaying NO LIMIT CYCLE
    text(-0.5, 0.5, 'NO LIMIT CYCLE', 'FontSize', 11, 'FontWeight', 'bold', ...
         'Color', 'r', 'BackgroundColor', [1 0.9 0.9], 'EdgeColor', 'r', 'Margin', 4);

    % Plot Legend
    legend('Aircraft Transfer Function: Gs\_ldynac(j\omega)', ...
           'Inverse Negative Describing Function: (-1 / N(A, j\omega))', ...
           'Critical Point (-1, 0i)', ...
           'Location', 'best');
end
hold off;

%% 6. Command Window Summary Output
fprintf('\n------------------- INTERSECTION PARAMETERS -----------------\n');
if HAS_INTERSECTION
    fprintf('\n=================== RATE LIMITING ELEMENT ===================\n');
    fprintf('  Intersection Point (x_int, y_int) : (%.2f, %.2f)\n', x_int, y_int);
    fprintf('  w_osc (Limit Cycle Frequency)    : %.2f rad/s\n', w_osc);
    fprintf('  A_rle                            : %.2f deg\n', A_rle);
    fprintf('\n============ AIRCRAFT LINEAR TRANSFER FUNCTION ==============\n');
    fprintf('  w_ac  (Plant Frequency)          : %.2f rad/s\n', w_osc);
    fprintf('  A_ac  (Plant Magnitude)          : %.2f\n', A_ac);
    fprintf('=============================================================\n');
else
    fprintf('  NO LIMIT CYCLE DETECTED (min_dist = %.4f >= tol = %.2f)\n', min_dist, tol);
    fprintf('=============================================================\n');
end