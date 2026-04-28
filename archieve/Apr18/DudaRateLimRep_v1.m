%% ========================================================================
%  DudaRateLimRep.m
%  Pitch-Attitude-Hold Autopilot with Rate Limiter Describing Function
%  ========================================================================
%  Reference: Duda — Example 4.6-1 (Transport Aircraft, 25000 ft, 500 ft/s)
%  Hanke (1994) AGARD-335 Rate Limiter Describing Function
%  ========================================================================

clear; clc; close all;

%% 1. Plant Matrices (Example 4.6-1)
ap = [-0.0082354,  18.938,      -32.170,     0.0,       5.9022e-05;
      -0.00025617, -0.56761,      0.0,        1.0,       2.2633e-06;
       0.0,         0.0,          0.0,        1.0,       0.0;
       1.3114e-05, -1.4847,       0.0,       -0.47599,  -1.4947e-07;
       0.0,        -500.00,     500.00,       0.0,       0.0];

bp = [0, 0, 0, -0.019781, 0]';
cp = [0, 0, 57.296, 0, 0;      % theta [deg]
      0, 0, 0,      57.296, 0]; % q     [deg/s]
dp = [0; 0];

plant = ss(ap, bp, cp, dp);

%% 2. Actuator (first-order lag, tau = 0.1 s, sign change)
aa = [-10];  ba = [10];
ca = [-1];   da = [0];
actua = ss(aa, ba, ca, da);

%% 3. Series Connection: Actuator -> Plant
sys1 = series(actua, plant);
[a, b, c, d] = ssdata(sys1);

%% 4. Autopilot Gains (Duda Example 4.6-1)
k_theta = 4.0;      % [deg elevator / deg pitch]
k_q     = 2.5;      % [deg elevator / (deg/s pitch rate)]

% Open-loop transfer function: L(s) = (k_theta + k_q*s) * G_theta(s)
%   G_theta = theta/delta_e = sys1(1,1)
%   G_q     = q/delta_e     = sys1(2,1)
%
% L(s) = k_theta * G_theta(s) + k_q * G_q(s)
G_theta = tf(sys1(1,1));
G_q     = tf(sys1(2,1));
L = k_theta * G_theta + k_q * G_q;    % Loop transfer function

fprintf('Open-Loop Transfer Function L(s) = (k_theta + k_q*s) * G_theta(s):\n');
L

% Closed-loop (for verification against Eq. 3 in textbook)
T_cl = feedback(L, 1);
fprintf('Closed-Loop Transfer Function theta/theta_c:\n');
T_cl

%% 5. Rate Limiter Parameters
R = 30;          % Rate limit [deg/s]
A = 10;          % Input amplitude [deg]
omega_onset = R / A;

fprintf('\n--- Rate Limiter Parameters ---\n');
fprintf('Rate limit R         = %.1f deg/s\n', R);
fprintf('Input amplitude A    = %.1f deg\n', A);
fprintf('Onset frequency      = %.2f rad/s\n', omega_onset);
fprintf('Full saturation at   = %.2f rad/s  (k=0.537)\n', omega_onset/0.537);

%% 6. Frequency Response Evaluation
w = logspace(-2, 1.5, 2000);

% --- L(jw): Linear open-loop ---
[mag_L, phase_L] = bode(L, w);
mag_L   = squeeze(mag_L);
phase_L = squeeze(phase_L);

% --- RL(A, w): Hanke describing function ---
k_vec = R ./ (A * w);
[gain_RL, phase_RL] = hanke_df(k_vec);

% --- Combined: RL * L ---
mag_comb   = gain_RL .* mag_L;
phase_comb = phase_RL + phase_L;

% --- dB ---
mag_L_dB    = 20*log10(mag_L);
mag_comb_dB = 20*log10(mag_comb);

% --- OLOP: at omega_onset ---
[~, idx_onset] = min(abs(w - omega_onset));
olop_phase = phase_comb(idx_onset);
olop_mag   = mag_comb_dB(idx_onset);

%% 7. BODE PLOT
figure('Name','Bode Plot','Position',[50 400 900 650]);

subplot(2,1,1);
semilogx(w, mag_L_dB, 'b-', 'LineWidth', 1.8); hold on;
semilogx(w, mag_comb_dB, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.8);
xline(omega_onset, '-', '\omega_{onset}', ...
    'Color', [0.1 0.6 0.3], 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
ylabel('Magnitude [dB]', 'FontSize', 11);
title(sprintf('Bode Plot:  k_\\theta = %.1f,  k_q = %.1f,  R = %.0f deg/s,  A = %.0f deg', ...
    k_theta, k_q, R, A), 'FontSize', 12);
legend('Linear L(s)', 'RL \times L(s)  (Hanke DF)', ...
    'Location', 'southwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10);

subplot(2,1,2);
semilogx(w, phase_L, 'b-', 'LineWidth', 1.8); hold on;
semilogx(w, phase_comb, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.8);
xline(omega_onset, '-', '', ...
    'Color', [0.1 0.6 0.3], 'LineWidth', 1.5);
ylabel('Phase [deg]', 'FontSize', 11);
xlabel('Frequency [rad/s]', 'FontSize', 11);
legend('Linear L(s)', 'RL \times L(s)  (Hanke DF)', ...
    'Location', 'southwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10);

%% 8. NYQUIST DIAGRAM
L_jw    = mag_L    .* exp(1j * phase_L    * pi/180);
comb_jw = mag_comb .* exp(1j * phase_comb * pi/180);

figure('Name','Nyquist Diagram','Position',[50 50 750 650]);
plot(real(L_jw), imag(L_jw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(comb_jw), imag(comb_jw), '--', 'Color', [0.8 0.1 0.1], ...
    'LineWidth', 1.8);
plot(-1, 0, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [0.6 0.0 0.6], ...
    'MarkerFaceColor', [0.7 0.2 0.7], 'LineWidth', 1.5);
% OLOP on Nyquist
olop_complex = mag_comb(idx_onset) * exp(1j * phase_comb(idx_onset) * pi/180);
plot(real(olop_complex), imag(olop_complex), 'o', ...
    'MarkerSize', 10, 'MarkerFaceColor', [0.1 0.6 0.3], ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
xlabel('Real Axis', 'FontSize', 11);
ylabel('Imaginary Axis', 'FontSize', 11);
title('Nyquist Diagram', 'FontSize', 12);
legend('Linear L(s)', 'RL \times L(s)', 'Critical Point (-1, 0)', ...
    sprintf('OLOP (%.1f rad/s)', omega_onset), ...
    'Location', 'best', 'FontSize', 10);
grid on; axis equal;
set(gca, 'FontSize', 10);

%% 9. NICHOLS CHART (with M-circles)
figure('Name','Nichols Chart','Position',[500 100 850 700]);
hold on;

% --- M-contours ---
M_dB_vals   = [-12 -6 -3 0 3 6];
M_labels    = {'-12 dB','-6 dB','-3 dB','0 dB','3 dB','6 dB'};
label_phases = [-120, -130, -140, -180, -200, -210];

for mi = 1:length(M_dB_vals)
    M_dB = M_dB_vals(mi);
    M    = 10^(M_dB/20);

    if abs(M - 1) < 1e-10
        plot([-180 -180], [-30 30], '-', ...
            'Color', [0.78 0.78 0.78], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
        text(-178, 13, '0 dB', 'FontSize', 9, 'Color', [0.55 0.55 0.55]);
        continue;
    end

    theta_mc = linspace(-359, -1, 4000) * pi/180;
    mag_mc = []; ph_mc = [];

    for ti = 1:length(theta_mc)
        phi = theta_mc(ti);
        aq = 1 - M^2;
        bq = -2 * M^2 * cos(phi);
        cq = -M^2;
        disc = bq^2 - 4*aq*cq;
        if disc < 0, continue; end
        x1 = (-bq + sqrt(disc)) / (2*aq);
        x2 = (-bq - sqrt(disc)) / (2*aq);
        for x = [x1 x2]
            if x > 1e-6
                mag_mc(end+1) = 20*log10(x); %#ok
                ph_mc(end+1)  = phi * 180/pi; %#ok
            end
        end
    end

    if ~isempty(mag_mc)
        plot(ph_mc, mag_mc, '-', 'Color', [0.78 0.78 0.78], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
        [~, li] = min(abs(ph_mc - label_phases(mi)));
        if ~isempty(li)
            text(ph_mc(li(1))+2, mag_mc(li(1))+0.5, M_labels{mi}, ...
                'FontSize', 9, 'Color', [0.55 0.55 0.55]);
        end
    end
end

% --- Linear L(s) ---
plot(phase_L, mag_L_dB, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'Linear');

% --- Describing Function RL x L(s) ---
plot(phase_comb, mag_comb_dB, '--', 'Color', [0.8 0.1 0.1], ...
    'LineWidth', 2, 'DisplayName', 'Describing Function');

% --- OLOP marker (green filled) ---
plot(olop_phase, olop_mag, 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0.0 0.0 0.0], ...
    'MarkerFaceColor', [0.1 0.6 0.3], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('OLOP  %.1f rad/s', omega_onset));

% --- Critical point (purple star) ---
plot(-180, 0, 'p', 'MarkerSize', 18, ...
    'MarkerEdgeColor', [0.5 0.0 0.5], ...
    'MarkerFaceColor', [0.7 0.2 0.7], 'LineWidth', 1.5, ...
    'DisplayName', 'Critical Point');

% --- Formatting ---
xlabel('Phase [deg]', 'FontSize', 12);
ylabel('Amplitude [dB]', 'FontSize', 12);
title(sprintf(['Nichols Chart: Open Loop Frequency Response\n' ...
    'k_\\theta = %.1f,  k_q = %.1f,  R = %.0f deg/s,  A = %.0f deg'], ...
    k_theta, k_q, R, A), 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
grid on;
xlim([-300 -50]);
ylim([-20 15]);
set(gca, 'XTick', -300:25:-50, 'FontSize', 10);

fprintf('\n--- OLOP ---\n');
fprintf('Phase = %.1f deg\n', olop_phase);
fprintf('Gain  = %.2f dB\n', olop_mag);
fprintf('at omega = %.2f rad/s\n', omega_onset);
fprintf('\nScript complete.\n');

%% === Local Function =====================================================
function [gain, phase_deg] = hanke_df(k)
    gain      = ones(size(k));
    phase_deg = zeros(size(k));

    % Region III (k <= 0.537): fully saturated
    idx3 = (k <= 0.537);
    gain(idx3)      = k(idx3) * pi/2;
    phase_deg(idx3) = -acosd( min(1, k(idx3)*pi/2) );

    % Region II (0.537 < k < 0.725): partially active
    idx2 = (k > 0.537) & (k < 0.725);
    gain(idx2)      = 1 - 4.51*(0.725 - k(idx2)).^2;
    phase_deg(idx2) = -173*(k(idx2) - 0.537) - 32.5;

    % Region I (k >= 0.725): inactive
end
