%% ========================================================================
%  DudaRateLimRep.m
%  Reproduction of OLOP Nichols Chart (Fig. 3 style)
%  ========================================================================
%  Plant:  Duda Example 4.6-1 (Transport, 25000 ft, 500 ft/s, xcg=0.25c)
%  RL:     Hanke (1994) AGARD-335 Describing Function
%
%  R = 34 deg/s, A = 10 deg selected to place OLOP at 0 dB crossing
%  for maximum visual similarity to the reference Figure 3.
%  ========================================================================

clear; clc; close all;

%% --- Color Definitions -------------------------------------------------
col_linear = [0.00 0.25 0.70];     % blue
col_df     = [0.80 0.10 0.10];     % red
col_olop   = [0.10 0.60 0.30];     % green
col_crit   = [0.55 0.00 0.55];     % purple
col_onset  = [0.90 0.55 0.00];     % amber

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

%% 3. Rate Limiter Parameters (optimized for Fig. 3 match)
R = 34;          % Rate limit [deg/s]
A = 10;          % Input amplitude [deg]
omega_onset = R / A;   % 3.4 rad/s  ->  OLOP lands at ~0 dB

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

%% 5. BODE PLOT
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

%% 6. NYQUIST DIAGRAM
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

%% 7. NICHOLS CHART — Figure 3 Reproduction
figure('Name','Nichols Chart','Position',[400 80 850 720],'Color','w');
ax = axes; hold(ax,'on');

% --- M-circle contours with shading ---
M_dB_list = [0.25 0.5 1 3 6 12];
shade_rgb = [0.91 0.91 0.91];

for mi = 1:length(M_dB_list)
    M_pos = 10^( M_dB_list(mi)/20);
    M_neg = 10^(-M_dB_list(mi)/20);
    for M = [M_pos, M_neg]
        [ph_c, mag_c] = nichols_contour(M);
        if M > 1
            fill(ph_c, mag_c, shade_rgb, ...
                'EdgeColor', [0.65 0.65 0.65], 'LineWidth', 0.5, ...
                'FaceAlpha', 0.35, 'HandleVisibility', 'off');
        else
            plot(ph_c, mag_c, '-', 'Color', [0.65 0.65 0.65], ...
                'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end
end

% 0 dB vertical
plot([-180 -180], [-40 40], '-', 'Color', [0.65 0.65 0.65], ...
    'LineWidth', 0.5, 'HandleVisibility', 'off');

% M-circle labels
text(-155,  0.5, '0 dB',  'FontSize',9, 'Color',[0.4 0.4 0.4]);
text(-150,  5,   '3 dB',  'FontSize',9, 'Color',[0.4 0.4 0.4]);
text(-160,  9,   '6 dB',  'FontSize',9, 'Color',[0.4 0.4 0.4]);
text(-110, -1.5, '-3 dB', 'FontSize',9, 'Color',[0.4 0.4 0.4]);

% --- Linear L(s) (blue, solid) ---
plot(phase_L, mag_L_dB, '-', 'Color', col_linear, ...
    'LineWidth', 2.2, 'DisplayName', 'Linear');

% --- Describing Function (red, dash-dot) ---
plot(phase_comb, mag_comb_dB, '-.', 'Color', col_df, ...
    'LineWidth', 2.2, 'DisplayName', 'Describing Function');

% --- OLOP marker (green filled circle) ---
plot(olop_phase, olop_mag, 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_olop, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Open-Loop Onset Point (OLOP)');

% --- Annotate OLOP ---
text(olop_phase + 5, olop_mag + 1.5, ...
    sprintf('%.1f rad/s', omega_onset), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', col_olop);

% --- Formatting ---
xlabel('Phase, deg', 'FontSize', 12);
ylabel('Amplitude, dB', 'FontSize', 12);
title(sprintf(['Nichols Chart: Open Loop Frequency Response\n' ...
    'R = %g deg/s,  A = %g deg,  \\omega_{onset} = %.1f rad/s'], ...
    R, A, omega_onset), 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
grid on;
xlim([-300 -50]);
ylim([-20 15]);
set(gca, 'XTick', -300:50:-50, 'FontSize', 10);
set(gca, 'Layer', 'top');

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
%  LOCAL FUNCTIONS
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

function [ph_deg, mag_dB] = nichols_contour(M)
    phi = linspace(-359.5, -0.5, 5000) * pi/180;
    ph_deg = [];  mag_dB = [];

    for i = 1:length(phi)
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

    [ph_deg, idx] = sort(ph_deg);
    mag_dB = mag_dB(idx);

    keep = (mag_dB > -30) & (mag_dB < 30);
    ph_deg = ph_deg(keep);
    mag_dB = mag_dB(keep);
end
