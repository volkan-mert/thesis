%% ========================================================================
%  DudaRateLimRep.m
%  Reproduction of OLOP Nichols Chart (Fig. 3 style)
%  ========================================================================
%  Plant:  Duda Example 4.6-1 (Transport, 25000 ft, 500 ft/s, xcg=0.25c)
%  RL:     Hanke (1994) AGARD-335 Describing Function
%  ========================================================================

clear; clc; close all;

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
L       = k_theta * G_theta + k_q * G_q;   % Open-loop TF

%% 3. Rate Limiter Parameters
R = 30;          % Rate limit [deg/s]
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

% OLOP index
[~, idx_onset] = min(abs(w - omega_onset));
olop_phase = phase_comb(idx_onset);
olop_mag   = mag_comb_dB(idx_onset);

%% 5. Bode Plot
figure('Name','Bode Plot','Position',[50 450 900 550], ...
    'Color','w');

subplot(2,1,1);
semilogx(w, mag_L_dB, '-', 'Color', [0.0 0.0 0.0], ...
    'LineWidth', 1.5); hold on;
semilogx(w, mag_comb_dB, '--', 'Color', [0.0 0.0 0.0], ...
    'LineWidth', 1.5);
xline(omega_onset, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
ylabel('Amplitude, dB', 'FontSize', 11);
title(sprintf('Bode Plot  (R = %g deg/s,  A = %g deg,  \\omega_{onset} = %g rad/s)', ...
    R, A, omega_onset), 'FontSize', 12);
legend('Linear', 'Describing Function', '\omega_{onset}', ...
    'Location','southwest', 'FontSize', 9);
grid on;  set(gca,'FontSize',10);

subplot(2,1,2);
semilogx(w, phase_L, '-', 'Color', [0.0 0.0 0.0], ...
    'LineWidth', 1.5); hold on;
semilogx(w, phase_comb, '--', 'Color', [0.0 0.0 0.0], ...
    'LineWidth', 1.5);
xline(omega_onset, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
ylabel('Phase, deg', 'FontSize', 11);
xlabel('Frequency, rad/s', 'FontSize', 11);
legend('Linear', 'Describing Function', ...
    'Location','southwest', 'FontSize', 9);
grid on;  set(gca,'FontSize',10);

%% 6. Nyquist Diagram
L_jw    = mag_L    .* exp(1j * phase_L    * pi/180);
comb_jw = mag_comb .* exp(1j * phase_comb * pi/180);

olop_cplx = mag_comb(idx_onset) * exp(1j*phase_comb(idx_onset)*pi/180);

figure('Name','Nyquist Diagram','Position',[50 50 700 600], ...
    'Color','w');
plot(real(L_jw), imag(L_jw), '-', 'Color', [0 0 0], ...
    'LineWidth', 1.5); hold on;
plot(real(comb_jw), imag(comb_jw), '--', 'Color', [0 0 0], ...
    'LineWidth', 1.5);
plot(-1, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
plot(real(olop_cplx), imag(olop_cplx), 'ko', ...
    'MarkerSize', 9, 'MarkerFaceColor', 'k');
xlabel('Real Axis', 'FontSize', 11);
ylabel('Imaginary Axis', 'FontSize', 11);
title('Nyquist Diagram', 'FontSize', 12);
legend('Linear', 'Describing Function', ...
    'Critical Point (-1, 0)', ...
    sprintf('OLOP (%.1f rad/s)', omega_onset), ...
    'Location','best', 'FontSize', 9);
grid on; axis equal;
set(gca,'FontSize',10);

%% 7. Nichols Chart — Figure 3 Reproduction
figure('Name','Nichols Chart','Position',[400 80 850 720], ...
    'Color','w');
ax = axes; hold(ax, 'on');

% ------------------------------------------------------------------
%  M-circle contours with shading (closed-loop magnitude)
% ------------------------------------------------------------------
M_dB_list = [0.25 0.5 1 3 6 12];   % positive dB values
shade_rgb = [0.91 0.91 0.91];       % light gray fill

% Draw symmetric M-circle pairs: +M dB and -M dB
for mi = 1:length(M_dB_list)
    M_pos = 10^( M_dB_list(mi)/20);
    M_neg = 10^(-M_dB_list(mi)/20);
    for M = [M_pos, M_neg]
        [ph_c, mag_c] = nichols_contour(M);
        if M > 1
            % Fill inside for positive-dB contours (upper lobe)
            fill(ph_c, mag_c, shade_rgb, ...
                'EdgeColor', [0.65 0.65 0.65], 'LineWidth', 0.5, ...
                'FaceAlpha', 0.35, 'HandleVisibility', 'off');
        else
            plot(ph_c, mag_c, '-', 'Color', [0.65 0.65 0.65], ...
                'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end
end

% 0 dB contour (vertical line at -180 deg)
plot([-180 -180], [-40 40], '-', 'Color', [0.65 0.65 0.65], ...
    'LineWidth', 0.5, 'HandleVisibility', 'off');

% M-circle labels
text(-155, 0.5,  '0 dB',  'FontSize', 9, 'Color', [0.4 0.4 0.4]);
text(-150, 5,    '3 dB',  'FontSize', 9, 'Color', [0.4 0.4 0.4]);
text(-160, 9,    '6 dB',  'FontSize', 9, 'Color', [0.4 0.4 0.4]);
text(-110, -1.5, '-3 dB', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

% ------------------------------------------------------------------
%  Linear open-loop L(s)
% ------------------------------------------------------------------
plot(phase_L, mag_L_dB, '-', 'Color', [0 0 0], ...
    'LineWidth', 2, 'DisplayName', 'Linear');

% ------------------------------------------------------------------
%  Describing Function  RL × L(s)
% ------------------------------------------------------------------
plot(phase_comb, mag_comb_dB, '-.', 'Color', [0 0 0], ...
    'LineWidth', 2, 'DisplayName', 'Describing Function');

% ------------------------------------------------------------------
%  OLOP marker
% ------------------------------------------------------------------
plot(olop_phase, olop_mag, 'ko', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Open-Loop Onset Point (OLOP)');

% Annotate OLOP with omega value
text(olop_phase + 5, olop_mag + 1.2, ...
    sprintf('%.1f rad/s', omega_onset), ...
    'FontSize', 10, 'FontWeight', 'bold');

% ------------------------------------------------------------------
%  Formatting
% ------------------------------------------------------------------
xlabel('Phase, deg', 'FontSize', 12);
ylabel('Amplitude, dB', 'FontSize', 12);
title('Nichols Chart: Open Loop Frequency Response', 'FontSize', 13);
legend('Location', 'southwest', 'FontSize', 10);
grid on;
xlim([-300 -50]);
ylim([-20 15]);
set(gca, 'XTick', -300:50:-50, 'FontSize', 10);
set(gca, 'Layer', 'top');

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
%  LOCAL FUNCTIONS
%  ========================================================================

function [gain, phase_deg] = hanke_df(k)
% HANKE_DF  Hanke (1994) rate limiter describing function.
%   k = R / (A*omega)
%   Region I   (k >= 0.725):  gain = 1,  phase = 0
%   Region II  (0.537 < k < 0.725):  polynomial fit
%   Region III (k <= 0.537):  fully saturated (peak ratio)

    gain      = ones(size(k));
    phase_deg = zeros(size(k));

    % Region III
    m3 = (k <= 0.537);
    gain(m3)      = k(m3) * pi/2;
    phase_deg(m3) = -acosd( min(1, k(m3)*pi/2) );

    % Region II
    m2 = (k > 0.537) & (k < 0.725);
    gain(m2)      = 1 - 4.51*(0.725 - k(m2)).^2;
    phase_deg(m2) = -173*(k(m2) - 0.537) - 32.5;
end

function [ph_deg, mag_dB] = nichols_contour(M)
% NICHOLS_CONTOUR  Compute a single closed-loop M-circle in Nichols coords.
%   M = linear magnitude of closed-loop |T(jw)|.
%   Returns phase (deg) and open-loop magnitude (dB) arrays.

    phi = linspace(-359.5, -0.5, 5000) * pi/180;
    ph_deg = [];
    mag_dB = [];

    for i = 1:length(phi)
        % Solve |G/(1+G)| = M  for |G|, given angle(G) = phi
        % |G|^2 (1 - M^2) - 2 M^2 cos(phi) |G| - M^2 = 0
        aq = 1 - M^2;
        bq = -2 * M^2 * cos(phi(i));
        cq = -M^2;
        disc = bq^2 - 4*aq*cq;
        if disc < 0, continue; end
        roots_x = [(-bq + sqrt(disc))/(2*aq), ...
                    (-bq - sqrt(disc))/(2*aq)];
        for x = roots_x
            if x > 1e-8
                ph_deg(end+1) = phi(i) * 180/pi; %#ok
                mag_dB(end+1) = 20*log10(x);      %#ok
            end
        end
    end

    % Sort for clean contour plotting
    [ph_deg, idx] = sort(ph_deg);
    mag_dB = mag_dB(idx);

    % Clip to visible window
    keep = (mag_dB > -30) & (mag_dB < 30);
    ph_deg = ph_deg(keep);
    mag_dB = mag_dB(keep);
end
