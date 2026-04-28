%% ========================================================================
%  DudaRateLimRep.m
%  Reproduction of OLOP Nichols Chart (Fig. 3 style)
%  ========================================================================
%  Plant:  Duda Example 4.6-1 (Transport, 25000 ft, 500 ft/s, xcg=0.25c)
%  RL:     Hanke (1994) AGARD-335 Describing Function
%  R = 34 deg/s, A = 10 deg  =>  OLOP at ~0 dB
%  ========================================================================

clear; clc; close all;

%% --- Color Definitions -------------------------------------------------
col_linear = [0.00 0.25 0.70];     % blue  — linear
col_df     = [0.80 0.10 0.10];     % red   — describing function
col_olop   = [0.10 0.60 0.30];     % green — OLOP
col_crit   = [0.00 0.30 0.85];     % blue  — critical point
col_onset  = [0.80 0.00 0.80];     % magenta — omega_onset

%% 1. Plant + Actuator
ap = [-0.0082354,  18.938,      -32.170,     0.0,       5.9022e-05;
      -0.00025617, -0.56761,      0.0,        1.0,       2.2633e-06;
       0.0,         0.0,          0.0,        1.0,       0.0;
       1.3114e-05, -1.4847,       0.0,       -0.47599,  -1.4947e-07;
       0.0,        -500.00,     500.00,       0.0,       0.0];

bp = [0, 0, 0, -0.019781, 0]';
cp = [0, 0, 57.296, 0, 0;
      0, 0, 0,      57.296, 0];
dp = [0; 0];

plant = ss(ap, bp, cp, dp);

aa = [-10]; ba = [10];
ca = [-1];  da = [0];
actua = ss(aa, ba, ca, da);

sys1 = series(actua, plant);
[a, b, c, d] = ssdata(sys1);

%% 2. Autopilot Gains
k_theta = 4.0;
k_q     = 2.5;

G_theta = tf(sys1(1,1));
G_q     = tf(sys1(2,1));
L       = k_theta * G_theta + k_q * G_q;

%% 3. Rate Limiter Parameters
R = 34;          % Rate limit [deg/s]
A = 10;          % Input amplitude [deg]
omega_onset = R / A;

%% 4. Frequency Response
w = logspace(-2, 1.5, 3000);

[mag_L, phase_L] = bode(L, w);
mag_L   = squeeze(mag_L);
phase_L = squeeze(phase_L);

k_vec = R ./ (A * w);
[gain_RL, phase_RL] = hanke_df(k_vec);

mag_comb   = gain_RL .* mag_L;
phase_comb = phase_RL + phase_L;

mag_L_dB    = 20*log10(mag_L);
mag_comb_dB = 20*log10(mag_comb);

[~, idx_onset] = min(abs(w - omega_onset));
olop_phase = phase_comb(idx_onset);
olop_mag   = mag_comb_dB(idx_onset);

%% ========================================================================
%  FIGURE 1: BODE PLOT
%  ========================================================================
figure('Name','Bode Plot','Position',[50 450 900 550],'Color','w');

subplot(2,1,1);
semilogx(w, mag_L_dB, '-', 'Color', col_linear, ...
    'LineWidth', 1.8); hold on;
semilogx(w, mag_comb_dB, '--', 'Color', col_df, ...
    'LineWidth', 1.8);
xline(omega_onset, '-', '\omega_{onset}', ...
    'Color', col_onset, 'LineWidth', 1.5, ...
    'LabelVerticalAlignment','bottom', 'FontSize',10);
ylabel('Amplitude, dB', 'FontSize', 11);
title(sprintf('Bode Plot  (R = %g deg/s,  A = %g deg,  \\omega_{onset} = %.1f rad/s)', ...
    R, A, omega_onset), 'FontSize', 12);
legend('Linear', 'Describing Function', '\omega_{onset}', ...
    'Location','southwest', 'FontSize', 9);
grid on; set(gca,'FontSize',10);

subplot(2,1,2);
semilogx(w, phase_L, '-', 'Color', col_linear, ...
    'LineWidth', 1.8); hold on;
semilogx(w, phase_comb, '--', 'Color', col_df, ...
    'LineWidth', 1.8);
xline(omega_onset, '-', '', 'Color', col_onset, 'LineWidth', 1.5);
ylabel('Phase, deg', 'FontSize', 11);
xlabel('Frequency, rad/s', 'FontSize', 11);
legend('Linear', 'Describing Function', ...
    'Location','southwest', 'FontSize', 9);
grid on; set(gca,'FontSize',10);

%% ========================================================================
%  FIGURE 2: NYQUIST DIAGRAM
%  ========================================================================
L_jw    = mag_L    .* exp(1j * phase_L    * pi/180);
comb_jw = mag_comb .* exp(1j * phase_comb * pi/180);
olop_cplx = mag_comb(idx_onset) * exp(1j*phase_comb(idx_onset)*pi/180);

figure('Name','Nyquist Diagram','Position',[50 50 700 600],'Color','w');
plot(real(L_jw), imag(L_jw), '-', 'Color', col_linear, ...
    'LineWidth', 1.8); hold on;
plot(real(comb_jw), imag(comb_jw), '--', 'Color', col_df, ...
    'LineWidth', 1.8);
plot(-1, 0, 'p', 'MarkerSize', 16, ...
    'MarkerEdgeColor', col_crit, 'MarkerFaceColor', col_crit);
plot(real(olop_cplx), imag(olop_cplx), 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_olop, ...
    'LineWidth', 1.5);
xlabel('Real Axis', 'FontSize', 11);
ylabel('Imaginary Axis', 'FontSize', 11);
title('Nyquist Diagram', 'FontSize', 12);
legend('Linear', 'Describing Function', ...
    'Critical Point (-1, 0)', ...
    sprintf('OLOP (%.1f rad/s)', omega_onset), ...
    'Location','best', 'FontSize', 9);
grid on; axis equal;
set(gca,'FontSize',10);

%% ========================================================================
%  FIGURE 3: NICHOLS CHART — Both curves
%  ========================================================================
figure('Name','Nichols Chart','Position',[400 80 850 720],'Color','w');

% Use MATLAB's built-in nichols grid for clean M-circles
nichols(L, w);  hold on;   % plots linear + draws ngrid automatically
% Clear the auto-plotted curve (we will re-draw with our colors)
delete(findobj(gca, 'Type', 'Line', 'Color', 'b'));
cla; hold on;

% Draw ngrid background (M and N circles)
ngrid;

% Linear L(s) — blue solid
plot(phase_L, mag_L_dB, '-', 'Color', col_linear, ...
    'LineWidth', 2.2, 'DisplayName', 'Linear');

% Describing Function — red dash-dot (single line object)
plot(phase_comb, mag_comb_dB, '-.', 'Color', col_df, ...
    'LineWidth', 2.2, 'DisplayName', 'Describing Function');

% OLOP — green
plot(olop_phase, olop_mag, 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_olop, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Open-Loop Onset Point (OLOP)');

text(olop_phase + 5, olop_mag + 1.5, ...
    sprintf('%.1f rad/s', omega_onset), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', col_olop);

xlabel('Phase, deg', 'FontSize', 12);
ylabel('Amplitude, dB', 'FontSize', 12);
title(sprintf(['Nichols Chart: Open Loop Frequency Response\n' ...
    'R = %g deg/s,  A = %g deg,  \\omega_{onset} = %.1f rad/s'], ...
    R, A, omega_onset), 'FontSize', 12);
legend('Linear', 'Describing Function', 'OLOP', ...
    'Location', 'southwest', 'FontSize', 10);
xlim([-300 -50]);
ylim([-20 15]);

%% ========================================================================
%  FIGURE 4: NICHOLS CHART — Describing Function Only
%  ========================================================================
figure('Name','Nichols Chart — DF Only','Position',[500 80 850 720],'Color','w');
hold on;

ngrid;

% Describing Function only — red dash-dot
plot(phase_comb, mag_comb_dB, '-.', 'Color', col_df, ...
    'LineWidth', 2.5, 'DisplayName', 'Describing Function');

% OLOP
plot(olop_phase, olop_mag, 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_olop, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Open-Loop Onset Point (OLOP)');

text(olop_phase + 5, olop_mag + 1.5, ...
    sprintf('%.1f rad/s', omega_onset), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', col_olop);

xlabel('Phase, deg', 'FontSize', 12);
ylabel('Amplitude, dB', 'FontSize', 12);
title(sprintf(['Nichols Chart: Describing Function RL \\times L(s)\n' ...
    'R = %g deg/s,  A = %g deg,  \\omega_{onset} = %.1f rad/s'], ...
    R, A, omega_onset), 'FontSize', 12);
legend('Describing Function', 'OLOP', ...
    'Location', 'southwest', 'FontSize', 10);
xlim([-300 -50]);
ylim([-20 15]);

%% Print Summary
fprintf('============================================\n');
fprintf('  OLOP Analysis — Best Fit for Fig. 3\n');
fprintf('============================================\n');
fprintf('  k_theta       = %.1f\n', k_theta);
fprintf('  k_q           = %.1f\n', k_q);
fprintf('  R             = %.1f deg/s\n', R);
fprintf('  A             = %.1f deg\n', A);
fprintf('  omega_onset   = %.2f rad/s\n', omega_onset);
fprintf('--------------------------------------------\n');
fprintf('  OLOP Phase    = %.1f deg\n', olop_phase);
fprintf('  OLOP Gain     = %.2f dB\n', olop_mag);
fprintf('============================================\n');

%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================

function [gain, phase_deg] = hanke_df(k)
    gain      = ones(size(k));
    phase_deg = zeros(size(k));

    m3 = (k <= 0.537);
    gain(m3)      = k(m3) * pi/2;
    phase_deg(m3) = -acosd( min(1, k(m3)*pi/2) );

    m2 = (k > 0.537) & (k < 0.725);
    gain(m2)      = 1 - 4.51*(0.725 - k(m2)).^2;
    phase_deg(m2) = -173*(k(m2) - 0.537) - 32.5;
end
