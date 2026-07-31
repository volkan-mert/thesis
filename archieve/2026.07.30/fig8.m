%% APPLIED NONLINEAR CONTROL, SLOTINE, 1991, Chapter 5, Page 184
% Limit cycle detection for frequency-dependent describing functions
% by sketching Nyquist Plot of G(jw) and -1 / N(A, w)

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
R = 15;                       % Slew rate / Max actuator rate (deg/s)
A = 0.1;                    % Input amplitude (deg) (Limit: A < 1.19226575 at R=60)
w = logspace(-1, 2, 1000);    % Frequency vector from 0.1 to 100 rad/s

% Calculate onset and critical saturation frequencies
w_onset = R./A;              % Saturation onset = 15 rad/s
w_crit  = 1.862 * w_onset;    % Fully developed saturation limit = 27.93 rad/s

mag = zeros(1, length(w));
phi = zeros(1, length(w));

for k = 1:length(w)
    alpha = w(k) / w_onset;   
    % if alpha < 1
    %     % Region I: No saturation (Eq. 15)
    %     mag(k) = 1; phi(k) = 0;
    %     disp('Region I')
    % elseif alpha < 1.862
    %     % Region II: Transition via cubic spline interpolation (Eq. 16)
    %     mag(k) = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
    %     phi(k) = 0.5280*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171; % rad
    %     disp('Region II')
    % else
    %     % Region III: Fully developed saturation (Eq. 17)
    %     wbar = 1 / alpha;    
    %     mag(k) = (4 * wbar) / pi;
    %     phi(k) = -acos(pi * wbar / 2);                                  % rad
    %     disp('Region III')
    % end

        % Region III: Fully developed saturation (Eq. 17)
        wbar = 1 / alpha;    
        mag(k) = (4 * wbar) / pi;
        phi(k) = -acos(pi * wbar / 2);                                  % rad
        disp('Region III')

end

%% 3. FRD Model Creation & Plotting
N_jw = mag .* exp(1j * phi);
sys_rle = frd(N_jw, w);
indf_rle = -1/sys_rle; 

figure;
nyquistplot(Gs_ac, indf_rle);
title('FREQUENCY-DEPENDENT DESCRIBING FUNCTION', '(NYQUIST PLOT)');
grid on;

%% 4. Mark the critical point (-1, 0)
hold on;
plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
text(-1.05, 0.2, ' (-1, 0i)', 'Color', 'r', 'FontWeight', 'bold');
ylim([-2, 2]);
xlim([-2, 2]);

%% 5. Geometric Intersection Search & Plotting

% Extract complex frequency response arrays
Gs_ac_vec    = squeeze(freqresp(Gs_ac, w));
indf_rle_vec = squeeze(freqresp(indf_rle, w));

% Simple nested loop to find closest point between the two curves
min_dist = inf;
best_i   = 1; % Index for Gs_ac_vec
best_j   = 1; % Index for indf_rle_vec

for i = 1:length(w)
    for j = 1:length(w)
        d = abs(Gs_ac_vec(i) - indf_rle_vec(j));
        if d < min_dist
            min_dist = d;
            best_i   = i;
            best_j   = j;
        end
    end
end

% Set distance tolerance to confirm actual geometric intersection
tol = 0.08; 
HAS_INTERSECTION = (min_dist < tol);

if HAS_INTERSECTION
    % Raw values at closest indices
    x_int_ac  = real(Gs_ac_vec(best_i));
    y_int_ac  = imag(Gs_ac_vec(best_i));
    w_ac      = w(best_i);

    x_int_rle = real(indf_rle_vec(best_j));
    y_int_rle = imag(indf_rle_vec(best_j));
    w_rle     = w(best_j);

    % Amplitudes
    A_rle     = abs(indf_rle_vec(best_j));
    A_ac      = abs(Gs_ac_vec(best_i));

    % Single rounded intersection point coordinates
    x_int = round((x_int_ac + x_int_rle) / 2, 2);
    y_int = round((y_int_ac + y_int_rle) / 2, 2);

    % Plot single red star intersection point
    plot(x_int, y_int, 'r*', 'MarkerSize', 12);

    % Text box with w_ac and A_rle
    text_str = sprintf('  \\omega_{ac} = %.2f rad/s\n  A_{rle} = %.2f deg', w_ac, A_rle);
    text(x_int + 0.1, y_int - 0.15, text_str, 'FontSize', 9);

    % Legend with Intersection Point
    legend('Linear Transfer Function of CLAW * Aircraft (\theta_c / \theta): (G(j\omega))', ...
           'Inverse Negative Describing Function: (-1 / N(A, j\omega))', ...
           'Critical Point (-1, 0i)', ...
           'Intersection Point', ...
           'Location', 'best');
else
    % Text box displaying NO LIMIT CYCLE
    text(-0.5, 0.5, 'NO LIMIT CYCLE', 'FontSize', 11, 'FontWeight', 'bold', ...
         'Color', 'r', 'BackgroundColor', [1 0.9 0.9], 'EdgeColor', 'r', 'Margin', 4);

    % Legend WITHOUT Intersection Point
    legend('Linear Transfer Function of CLAW * Aircraft (\theta_c / \theta): (G(j\omega))', ...
           'Inverse Negative Describing Function: (-1 / N(A, j\omega))', ...
           'Critical Point (-1, 0i)', ...
           'Location', 'best');
end
hold off;

%% 6. Print Parameters in Command Window
fprintf('\n------------------- INTERSECTION PARAMETERS -----------------\n');
if HAS_INTERSECTION
    fprintf('\n=================== RATE LIMITING ELEMENT ===================\n');
    fprintf('  Intersection Point (x_int, y_int) : (%.2f, %.2f)\n', x_int, y_int);
    fprintf('  w_rle : %.2f rad/s\n', w_rle);
    fprintf('  A_rle : %.2f deg\n', A_rle);
    fprintf('\n============ AIRCRAFT LINEAR TRANSFER FUNCTION ==============\n');
    fprintf('  w_ac  : %.2f rad/s\n', w_ac);
    fprintf('  A_ac  : %.2f deg\n', A_ac);
    fprintf('=============================================================\n');
else
    fprintf('  NO LIMIT CYCLE DETECTED\n');
    fprintf('=============================================================\n');
end