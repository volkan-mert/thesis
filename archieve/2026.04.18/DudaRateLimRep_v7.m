%% ========================================================================
%  DudaRateLimRep.m
%  RL, G(s), and RL×G(s) — Bode, Nyquist, Nichols
%  ========================================================================
%  Plant:  Duda Example 4.6-1 (Transport, 25000 ft, 500 ft/s, xcg=0.25c)
%  RL:     Hanke (1994) AGARD-335 Describing Function
%  R = 34 deg/s, A = 10 deg  =>  omega_onset = 3.4 rad/s
%  ========================================================================

clear; clc; close all;

%% --- Color Definitions -------------------------------------------------
col_plant  = [0.00 0.25 0.70];     % blue    — G(s) linear
col_rl     = [0.80 0.00 0.80];     % magenta — RL only
col_rlg    = [0.80 0.10 0.10];     % red     — RL × G
col_olop   = [0.10 0.60 0.30];     % green   — OLOP
col_onset  = [0.90 0.55 0.00];     % amber   — omega_onset

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

%% 2. Autopilot Gains
k_theta = 4.0;
k_q     = 2.5;

G_theta = tf(sys1(1,1));
G_q     = tf(sys1(2,1));
L       = k_theta * G_theta + k_q * G_q;   % Open-loop with controller

%% 3. Rate Limiter Parameters
R = 34;
A = 10;
omega_onset = R / A;

%% 4. Frequency Response
w = logspace(-2, 1.5, 3000);

% --- G(s): plant + actuator + autopilot ---
[mag_G, phase_G] = bode(L, w);
mag_G   = squeeze(mag_G);
phase_G = squeeze(phase_G);
mag_G_dB = 20*log10(mag_G);

% --- RL: Hanke describing function ---
k_vec = R ./ (A * w);
[gain_RL, phase_RL] = hanke_df(k_vec);
gain_RL_dB = 20*log10(gain_RL);

% --- RL × G: combined ---
mag_comb   = gain_RL .* mag_G;
phase_comb = phase_RL + phase_G;
mag_comb_dB = 20*log10(mag_comb);

% --- OLOP ---
[~, idx_onset] = min(abs(w - omega_onset));
olop_phase = phase_comb(idx_onset);
olop_mag   = mag_comb_dB(idx_onset);

%% ========================================================================
%  FIGURE 1: BODE PLOT
%  ========================================================================
figure('Name','Bode Plot','Position',[50 450 900 600],'Color','w');

subplot(2,1,1);
semilogx(w, mag_G_dB, '-', 'Color', col_plant, ...
    'LineWidth', 1.8); hold on;
semilogx(w, gain_RL_dB, '-', 'Color', col_rl, ...
    'LineWidth', 1.8);
semilogx(w, mag_comb_dB, '--', 'Color', col_rlg, ...
    'LineWidth', 1.8);
xline(omega_onset, ':', 'Color', col_onset, 'LineWidth', 1.5);
ylabel('Magnitude, dB', 'FontSize', 11);
title(sprintf('Bode Plot  (R = %g deg/s,  A = %g deg,  \\omega_{onset} = %.1f rad/s)', ...
    R, A, omega_onset), 'FontSize', 12);
legend('G(s)  Linear', 'RL  Describing Function', ...
    'RL \times G(s)', '\omega_{onset}', ...
    'Location','southwest', 'FontSize', 9);
grid on; set(gca,'FontSize',10);

subplot(2,1,2);
semilogx(w, phase_G, '-', 'Color', col_plant, ...
    'LineWidth', 1.8); hold on;
semilogx(w, phase_RL, '-', 'Color', col_rl, ...
    'LineWidth', 1.8);
semilogx(w, phase_comb, '--', 'Color', col_rlg, ...
    'LineWidth', 1.8);
xline(omega_onset, ':', 'Color', col_onset, 'LineWidth', 1.5);
ylabel('Phase, deg', 'FontSize', 11);
xlabel('Frequency, rad/s', 'FontSize', 11);
legend('G(s)  Linear', 'RL  Describing Function', ...
    'RL \times G(s)', ...
    'Location','southwest', 'FontSize', 9);
grid on; set(gca,'FontSize',10);

%% ========================================================================
%  FIGURE 2: NYQUIST DIAGRAM
%  ========================================================================
G_jw    = mag_G    .* exp(1j * phase_G    * pi/180);
RL_jw   = gain_RL  .* exp(1j * phase_RL   * pi/180);
comb_jw = mag_comb .* exp(1j * phase_comb * pi/180);
olop_cplx = comb_jw(idx_onset);

figure('Name','Nyquist Diagram','Position',[50 50 750 650],'Color','w');
plot(real(G_jw), imag(G_jw), '-', 'Color', col_plant, ...
    'LineWidth', 1.8); hold on;
plot(real(RL_jw), imag(RL_jw), '-', 'Color', col_rl, ...
    'LineWidth', 1.8);
plot(real(comb_jw), imag(comb_jw), '--', 'Color', col_rlg, ...
    'LineWidth', 1.8);
plot(-1, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
plot(real(olop_cplx), imag(olop_cplx), 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_olop, ...
    'LineWidth', 1.5);
xlabel('Real Axis', 'FontSize', 11);
ylabel('Imaginary Axis', 'FontSize', 11);
title('Nyquist Diagram', 'FontSize', 12);
legend('G(s)  Linear', 'RL  Describing Function', ...
    'RL \times G(s)', 'Critical Point (-1, 0)', ...
    sprintf('OLOP (%.1f rad/s)', omega_onset), ...
    'Location','best', 'FontSize', 9);
grid on; axis equal;
set(gca,'FontSize',10);

%% ========================================================================
%  FIGURE 3: NICHOLS CHART
%  ========================================================================
figure('Name','Nichols Chart','Position',[400 80 850 720],'Color','w');
hold on;

ngrid;

% G(s) linear — blue solid
plot(phase_G, mag_G_dB, '-', 'Color', col_plant, ...
    'LineWidth', 2.2, 'DisplayName', 'G(s)  Linear');

% RL only — magenta solid
plot(phase_RL, gain_RL_dB, '-', 'Color', col_rl, ...
    'LineWidth', 2.2, 'DisplayName', 'RL  Describing Function');

% RL × G — red dash-dot
plot(phase_comb, mag_comb_dB, '-.', 'Color', col_rlg, ...
    'LineWidth', 2.2, 'DisplayName', 'RL \times G(s)');

% OLOP — green
plot(olop_phase, olop_mag, 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_olop, ...
    'LineWidth', 1.5, ...
    'DisplayName', sprintf('OLOP (%.1f rad/s)', omega_onset));

text(olop_phase + 5, olop_mag + 1.5, ...
    sprintf('%.1f rad/s', omega_onset), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', col_olop);

xlabel('Phase, deg', 'FontSize', 12);
ylabel('Amplitude, dB', 'FontSize', 12);
title(sprintf(['Nichols Chart: Open Loop Frequency Response\n' ...
    'R = %g deg/s,  A = %g deg,  \\omega_{onset} = %.1f rad/s'], ...
    R, A, omega_onset), 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
xlim([-300 -50]);
ylim([-25 15]);

%% Print Summary
fprintf('============================================\n');
fprintf('  OLOP Analysis Summary\n');
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
