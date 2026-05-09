clear;clc;close all

%% Example : Draw Fig. 16 Nichols chart: open-loop frequency response q∕qc 
% 2025, Rodrigues et al, CEAS Aeronautical Journal (2026), 
% "Conditional integrator sliding mode control to reduce susceptibility 
% to pilot‑induced oscillations"

num_claw_7 = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7 = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7 = tf(num_claw_7, den_claw_7);

num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7 = tf(num_ldynac_7, den_ldynac_7);

Gs_ac_7 = Gs_claw_7*Gs_ldynac_7; % TF of the open-loop by q / q_c.

% Picking the cursor coordinates on the graph

pick_cursor_magnitude_olop = 2.0629; % omega_OLOP Magnitude in dB
pick_cursor_phase_olop = -136.6216; % omega_OLOP Phase in degrees

%% Calculating \omega_{onset}

F_uc_2_urle = feedback(Gs_claw_7,Gs_ldynac_7); % TF of the open-loop by q / q_c.

% Double check the value of the closed-loop onset frequency
R = deg2rad(60);
w=5.077;
Ain = 13.68/180*pi;  % since pilot input is 100% and kf=13.68
rate = Ain*abs(evalfr(F_uc_2_urle,1j*w))*w;

% Find the closed-loop onset frequency
f = @(w) (R - Ain*abs(evalfr(F_uc_2_urle,1j*w))*w)^2;
A = []; b = []; Aeq = []; beq = [];
lb = 0; ub = 100;
w0 = 0;
[w_opt, fval] = fmincon(f, w0, A, b, Aeq, beq, lb, ub)

%% RLE_DescribingFunction.m
%  Piecewise describing function (DF) of a Rate Limiter Element (RLE).
%
%  Three-region formulation as cited in the paper (Eqs. 15-17),
%  attributed to Duda / Hanke (AGARD-AR-335) / Gilbreath:
%
%     ϖ      = ω_onset / ω
%     α      = ω / ω_onset                    (= ϖ^-1)
%     ω_onset= R / A_{i,RLE}                  (R is the rate limit)
%
%  Region I    ω/ω_onset <  1.000  : no saturation
%                                     A   = 1
%                                     φ   = 0
%
%  Region II   1 ≤ ω/ω_onset < 1.862 : transition (cubic spline in α)
%                                     A   = 0.2908α^3 - 1.4396α^2 + 1.9232α + 0.2230
%                                     φ   = 0.5280α^3 - 2.6213α^2 + 3.5056α - 1.4171   [rad]
%
%  Region III  ω/ω_onset ≥ 1.862  : fully developed saturation
%                                     A   = 4ϖ/π
%                                     φ   = -acos(πϖ/2)                                 [rad]
%
%  Outputs:
%     - Bode-style plot (matches Fig. 12 of the reference)
%     - Nichols-style plot of N(jω,A) alone — the locus whose shape is
%       reproduced by the red curve in Fig. 16 once it is combined with
%       the linear open-loop response L(jω).
%

%% --------- Parameters ---------------------------------------------------
omega_onset = 1;                           % normalized onset frequency [rad/s]
omega       = logspace(-1, 2, 4000);       % frequency sweep            [rad/s]
xi          = omega / omega_onset;         % ω/ω_onset

A   = zeros(size(omega));                  % |N(jω,A)|
phi = zeros(size(omega));                  % ∠N(jω,A) in radians

%% --------- Region I  (no saturation) ------------------------------------
m1        = xi < 1;
A(m1)     = 1;
phi(m1)   = 0;

%% --------- Region II (transition, cubic spline in α = ω/ω_onset) --------
m2        = (xi >= 1) & (xi < 1.862);
a         = xi(m2);
A(m2)     = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.2230;
phi(m2)   = 0.5280*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;

%% --------- Region III (fully saturated, ϖ = ω_onset/ω) ------------------
m3        = xi >= 1.862;
v         = omega_onset ./ omega(m3);
A(m3)     = 4*v/pi;
phi(m3)   = -acos(pi*v/2);

%% --------- Convert to Bode-friendly units -------------------------------
A_dB      = 20*log10(A);
phi_deg   = rad2deg(phi);

%% --------- Continuity check at the boundaries ---------------------------
[~, kI ] = min(abs(xi - 1.000));
[~, kII] = min(abs(xi - 1.862));
fprintf('\nBoundary check:\n');
fprintf('  ω/ω_onset = 1.000 :  |N| = %+6.3f dB,  ∠N = %+7.3f deg\n', ...
    A_dB(kI),  phi_deg(kI));
fprintf('  ω/ω_onset = 1.862 :  |N| = %+6.3f dB,  ∠N = %+7.3f deg\n\n', ...
    A_dB(kII), phi_deg(kII));

%% --------- Plot 1 : Bode of the describing function (≈ Fig. 12) ---------
figure('Color','w','Name','RLE DF — Bode','Position',[80 80 920 620]);

subplot(2,1,1);
semilogx(xi, A_dB, 'b', 'LineWidth', 1.8); hold on; grid on; box on;
xline(1.000, '--k', 'LineWidth', 0.9);
xline(1.862, '--k', 'LineWidth', 0.9);
ylabel('|N(j\omega,A)|  [dB]');
title('RLE describing function — three-region formulation');
text(1.02,    2, '\omega = \omega_{onset}',       'Rotation',90,'FontSize',9);
text(1.90,    2, '\omega = 1.862\,\omega_{onset}','Rotation',90,'FontSize',9);

subplot(2,1,2);
semilogx(xi, phi_deg, 'r', 'LineWidth', 1.8); hold on; grid on; box on;
xline(1.000, '--k', 'LineWidth', 0.9);
xline(1.862, '--k', 'LineWidth', 0.9);
xlabel('\omega/\omega_{onset}');
ylabel('\angle N(j\omega,A)  [deg]');
ylim([-100 10]);

%% --------- Plot 2 : Nichols-style locus of N(jω,A) ----------------------
figure('Color','w','Name','RLE DF — Nichols-style','Position',[80 80 760 600]);

plot(phi_deg(m1), A_dB(m1), 'Color',[0.20 0.65 0.30], 'LineWidth', 2.2); hold on;
plot(phi_deg(m2), A_dB(m2), 'Color',[0.58 0.40 0.74], 'LineWidth', 2.2);
plot(phi_deg(m3), A_dB(m3), 'Color',[0.84 0.15 0.16], 'LineWidth', 2.2);
grid on; box on;

% Boundary markers
plot(0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
text(-3, 1.6, '\omega = \omega_{onset}', 'FontSize', 10);

phi_b = rad2deg(-acos(pi*0.537/2));
A_b   = 20*log10(4*0.537/pi);
plot(phi_b, A_b, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
text(phi_b+1, A_b+1.2, '\omega = 1.862\,\omega_{onset}', 'FontSize', 10);

xlabel('Phase  \angle N(j\omega,A)  [deg]');
ylabel('Gain  |N(j\omega,A)|  [dB]');
title('RLE describing function on Nichols-style axes');
legend({'Region I  (no saturation)', ...
        'Region II (transition, cubic)', ...
        'Region III (fully saturated)'}, ...
        'Location','southwest');
% xlim([-95 5]);
% ylim([-40 5]);

% Holding on Nichols Chart
hold on
%
opts = nicholsoptions;
opts.PhaseMatching      = 'on';
opts.PhaseMatchingFreq  = 1;        % rad/s
opts.PhaseMatchingValue = -180;     % anchor phase at that frequency
opts.PhaseWrapping = 'on';
opts.PhaseWrappingBranch = -360; % Wraps between -360 and 0
nicholsplot(Gs_ac_7, opts); 
hold on
%
%
plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop,'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
hold on
text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, sprintf('\\omega_{OLOP} = %.3f', w_opt), 'Color', 'red', 'FontWeight', 'bold', 'FontSize', 12);
%
grid on
% xlim([-300 -50]);
% ylim([-20 15]);
xlim([-300 5]);
ylim([-40 15]);

%% --------- Optional : combine with a linear plant to reproduce Fig. 16 --
% Uncomment and provide your closed-loop F_qc^δc(jω) (or any L(jω))
% to see the red curve of Fig. 16 emerge:
%
%   s    = 1j*omega;
%   L    = <your closed-loop transfer function evaluated on s>;
%   N    = A .* exp(1j*phi);
%   LN   = L .* N;
%   figure('Color','w');
%   plot(rad2deg(angle(LN)), 20*log10(abs(LN)), 'r', 'LineWidth', 1.6);
%   grid on; xlabel('Open-loop phase [deg]'); ylabel('Open-loop gain [dB]');
%   title('Open-loop response with RLE describing function (cf. Fig. 16)');
