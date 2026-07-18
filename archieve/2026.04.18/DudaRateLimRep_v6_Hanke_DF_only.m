%% ========================================================================
%  DudaRateLimRep.m
%  Hanke (1994) Rate Limiter Describing Function — RL Only
%  ========================================================================
%  R = 34 deg/s, A = 10 deg  =>  omega_onset = 3.4 rad/s
%  ========================================================================

clear; clc; close all;

%% --- Color Definitions -------------------------------------------------
col_df     = [0.80 0.10 0.10];     % red   — describing function
col_olop   = [0.10 0.60 0.30];     % green — OLOP
col_onset  = [0.80 0.00 0.80];     % magenta — omega_onset

%% 1. Rate Limiter Parameters
R = 34;          % Rate limit [deg/s]
A = 10;          % Input amplitude [deg]
omega_onset = R / A;

%% 2. Evaluate Hanke Describing Function
w = logspace(-1, 1.5, 3000);

k_vec = R ./ (A * w);
[gain_RL, phase_RL] = hanke_df(k_vec);

gain_RL_dB = 20*log10(gain_RL);

% OLOP index (at onset frequency)
[~, idx_onset] = min(abs(w - omega_onset));

%% ========================================================================
%  FIGURE 1: BODE PLOT — RL Only
%  ========================================================================
figure('Name','Bode Plot — RL','Position',[50 450 900 550],'Color','w');

subplot(2,1,1);
semilogx(w, gain_RL_dB, '-', 'Color', col_df, 'LineWidth', 2); hold on;
xline(omega_onset, '-', '\omega_{onset}', ...
    'Color', col_onset, 'LineWidth', 1.5, ...
    'LabelVerticalAlignment','bottom', 'FontSize',10);
ylabel('Gain, dB', 'FontSize', 11);
title(sprintf('Hanke Describing Function  (R = %g deg/s,  A = %g deg)', R, A), ...
    'FontSize', 12);
legend('RL(j\omega)', '\omega_{onset}', ...
    'Location','southwest', 'FontSize', 10);
grid on; set(gca,'FontSize',10);
% ylim([-25 2]);
ylim([-720 720]);

subplot(2,1,2);
semilogx(w, phase_RL, '-', 'Color', col_df, 'LineWidth', 2); hold on;
xline(omega_onset, '-', '', 'Color', col_onset, 'LineWidth', 1.5);
ylabel('Phase, deg', 'FontSize', 11);
xlabel('Frequency, rad/s', 'FontSize', 11);
legend('RL(j\omega)', 'Location','southwest', 'FontSize', 10);
grid on; set(gca,'FontSize',10);
% ylim([-100 5]);
ylim([-50 50]);

%% ========================================================================
%  FIGURE 2: NYQUIST DIAGRAM — RL Only
%  ========================================================================
RL_jw = gain_RL .* exp(1j * phase_RL * pi/180);

figure('Name','Nyquist — RL','Position',[50 50 700 600],'Color','w');
plot(real(RL_jw), imag(RL_jw), '-', 'Color', col_df, 'LineWidth', 2); hold on;
plot(real(RL_jw(idx_onset)), imag(RL_jw(idx_onset)), 'o', ...
    'MarkerSize', 10, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', col_olop, 'LineWidth', 1.5);
plot(1, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2);   % unity point
xlabel('Real', 'FontSize', 11);
ylabel('Imaginary', 'FontSize', 11);
title('Nyquist — Rate Limiter Describing Function', 'FontSize', 12);
legend('RL(j\omega)', sprintf('Onset (%.1f rad/s)', omega_onset), ...
    'Unity', 'Location','best', 'FontSize', 10);
grid on; axis equal;
set(gca,'FontSize',10);

%% ========================================================================
%  FIGURE 3: NICHOLS CHART — RL Only
%  ========================================================================
figure('Name','Nichols — RL','Position',[400 80 850 720],'Color','w');
hold on;

ngrid;

% RL curve
plot(phase_RL, gain_RL_dB, '-', 'Color', col_df, ...
    'LineWidth', 2.5, 'DisplayName', 'RL(j\omega)');

% Onset point
plot(phase_RL(idx_onset), gain_RL_dB(idx_onset), 'o', ...
    'MarkerSize', 12, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', col_olop, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Onset (%.1f rad/s)', omega_onset));

text(phase_RL(idx_onset) + 3, gain_RL_dB(idx_onset) + 1.2, ...
    sprintf('%.1f rad/s', omega_onset), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', col_olop);

xlabel('Phase, deg', 'FontSize', 12);
ylabel('Gain, dB', 'FontSize', 12);
title(sprintf(['Nichols Chart: Rate Limiter Describing Function\n' ...
    'R = %g deg/s,  A = %g deg'], R, A), 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
% xlim([-100 5]);
% ylim([-25 2]);
xlim([-720 720]);
ylim([-50 50]);

grid on

%% Print Summary
fprintf('============================================\n');
fprintf('  Hanke Describing Function — RL Only\n');
fprintf('============================================\n');
fprintf('  R             = %.1f deg/s\n', R);
fprintf('  A             = %.1f deg\n', A);
fprintf('  omega_onset   = %.2f rad/s\n', omega_onset);
fprintf('--------------------------------------------\n');
fprintf('  At onset:  Gain = %.2f dB,  Phase = %.1f deg\n', ...
    gain_RL_dB(idx_onset), phase_RL(idx_onset));
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
