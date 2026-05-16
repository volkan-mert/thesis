%% Fig16_OLOP_Nichols.m
% Reproduce Fig. 16 (Klyde / Duda OLOP-style Nichols chart):
%   blue  : linear open-loop response  CLAW * AC
%   red   : open-loop with rate-limiter describing function  CLAW * N_rle * AC
%   green : open-loop onset point (OLOP) at omega = omega_onset
%
% N_rle uses Duda's piecewise (cubic-spline transition) DF.
%
% Reference: Duda 1995 (AIAA-95-3304), AGARD-AR-335; Klyde-Mitchell.

clear; clc; close all;

%% Plant transfer functions
num_claw_7 = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7 = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7  = tf(num_claw_7, den_claw_7);

num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7  = tf(num_ldynac_7, den_ldynac_7);

Gs_lin = Gs_claw_7 * Gs_ldynac_7;        % CLAW*AC (linear OL)

%% Rate-limiter parameters
% Sweep A_i (or R) to move the red curve and the OLOP relative to the
% Nichols stability boundary. The reference Fig. 16 used omega_onset ~ 5 rad/s.
R   = 5.0;       % rate limit [u_rle / s]
A_i = 1.0;       % sinusoidal input amplitude into the RL
w_onset = R / A_i;

%% Frequency grid (densify across the three regions of the DF)
w_log  = logspace(-1, 2, 1000);
w_dens = linspace(0.6*w_onset, 3.0*w_onset, 80);
w = sort(unique([w_log, w_dens]));

%% Linear and DF open-loop responses
H_lin = squeeze(freqresp(Gs_lin, w));
N_rle = duda_rl_df(w, w_onset);
H_df  = H_lin(:) .* N_rle(:);

% Nichols coordinates
mag_lin_dB = 20*log10(abs(H_lin));
ph_lin_deg = unwrap(angle(H_lin)) * 180/pi;

mag_df_dB  = 20*log10(abs(H_df));
ph_df_deg  = unwrap(angle(H_df)) * 180/pi;

%% Open-loop onset point (OLOP)
% At w = w_onset the DF is exactly at the Region I/II boundary
% (A = 1, phi = 0), so the OLOP coincides with the linear response there.
H_olop      = squeeze(freqresp(Gs_lin, w_onset));
olop_mag_dB = 20*log10(abs(H_olop));
olop_ph_deg = unwrap([angle(H_lin(1)); angle(H_olop)]) * 180/pi;
olop_ph_deg = olop_ph_deg(end);

%% Plot — Nichols
figure('Color','w','Position',[80 80 900 680]);

% Nichols stability grid first (so curves draw on top)
ngrid;
hold on;

% Linear OL — blue with open circles (every Nth point to keep readable)
step_b = max(1, round(numel(w)/40));
plot(ph_lin_deg, mag_lin_dB, '-', 'Color',[0.10 0.35 0.85], 'LineWidth', 1.4, ...
     'DisplayName', 'Linear frequency response');
plot(ph_lin_deg(1:step_b:end), mag_lin_dB(1:step_b:end), 'o', ...
     'Color',[0.10 0.35 0.85], 'MarkerSize', 5, 'LineWidth', 1.0, ...
     'HandleVisibility','off');

% DF OL — red with filled dots
step_r = max(1, round(numel(w)/50));
plot(ph_df_deg, mag_df_dB, '-', 'Color',[0.85 0.10 0.10], 'LineWidth', 1.6, ...
     'DisplayName', 'Describing function');
plot(ph_df_deg(1:step_r:end), mag_df_dB(1:step_r:end), '.', ...
     'Color',[0.85 0.10 0.10], 'MarkerSize', 14, 'HandleVisibility','off');

% OLOP — filled green circle
plot(olop_ph_deg, olop_mag_dB, 'o', 'MarkerSize', 10, 'LineWidth', 1.2, ...
     'MarkerEdgeColor',[0 0.45 0], 'MarkerFaceColor',[0.20 0.85 0.20], ...
     'DisplayName','Open-loop onset point');

xlabel('Open-loop phase (^\circ)','FontSize',11);
ylabel('Open-loop gain (dB)','FontSize',11);
title(sprintf(['Open-loop CLAW \\cdot N_{rle} \\cdot AC   |   ' ...
               'R = %.3g, A_i = %.3g,  \\omega_{onset} = %.3f rad/s'], ...
               R, A_i, w_onset), 'FontSize', 11);
legend('Location','southwest','FontSize',10);
xlim([-300 -60]); ylim([-20 15]);
grid on; box on;

% Annotation showing the onset frequency near the green dot
text(olop_ph_deg+3, olop_mag_dB+1.2, ...
     sprintf('\\omega_{onset} = %.3f rad/s', w_onset), ...
     'Color',[0 0.45 0],'FontSize',9);

fprintf('\nOLOP coordinates:\n');
fprintf('  phase  = %7.2f deg\n', olop_ph_deg);
fprintf('  magn   = %7.2f dB\n',  olop_mag_dB);


%% ===== local function: Duda piecewise rate-limiter DF =====
function [N, A, phi] = duda_rl_df(w, w_onset)
% Returns complex DF N(jw, A_i) on the frequency vector w, given w_onset = R/A_i.
% Region I  : alpha < 1                     -> A=1,           phi=0
% Region II : 1 <= alpha < 1.862             -> cubic spline in alpha
% Region III: alpha >= 1.862                 -> A=4*wbar/pi,  phi=-acos(pi*wbar/2)
% with alpha = w/w_onset, wbar = 1/alpha.

    alpha = w(:).' / w_onset;
    A     = zeros(size(alpha));
    phi   = zeros(size(alpha));

    i1 = alpha < 1;
    A(i1)   = 1;
    phi(i1) = 0;

    i2 = (alpha >= 1) & (alpha < 1.862);
    a  = alpha(i2);
    A(i2)   = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.2230;
    phi(i2) = 0.5280*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;

    i3 = alpha >= 1.862;
    wbar = 1 ./ alpha(i3);
    A(i3)   = 4*wbar / pi;
    phi(i3) = -acos(pi*wbar/2);

    N = A .* exp(1j*phi);
end
