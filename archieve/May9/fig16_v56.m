clear; clc; close all;

%% =========================================================================
%  Reproduce Fig. 16 — Open-loop Nichols of q/q_c with RLE describing fn.
%  Reference: Rodrigues et al., CEAS Aeronautical Journal (2026)
%% =========================================================================

%% --------- Linear plant : control law × aircraft dynamics ----------------
num_claw_7   = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7    = tf(num_claw_7, den_claw_7);

num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7  = tf(num_ldynac_7, den_ldynac_7);

Gs_ac_7      = Gs_claw_7 * Gs_ldynac_7;          % open-loop  q / q_c

% Cursor-picked OLOP point from the published Fig. 16 (matched-phase axis)
pick_cursor_magnitude_olop = 2.0629;             % dB
pick_cursor_phase_olop     = -136.6216;          % deg

%% --------- Closed-loop onset frequency  ω_onset  (Eq. 18) ----------------
F_uc_2_urle = feedback(Gs_claw_7, Gs_ldynac_7);

R   = deg2rad(60);            % rate limit                 [rad/s]
Ain = 13.68/180*pi;           % pilot input × K_f          [rad]

f       = @(w) (R - Ain*abs(evalfr(F_uc_2_urle, 1j*w))*w)^2;
fminopt = optimoptions('fmincon','Display','off');
[w_opt, ~] = fmincon(f, 5, [],[],[],[], 0, 100, [], fminopt);

fprintf('Closed-loop onset frequency  ω_onset = %.4f rad/s\n\n', w_opt);

%% =========================================================================
%  Piecewise RLE describing function   (Duda / Hanke / Gilbreath)
%% =========================================================================
omega_onset = w_opt;
omega       = logspace(-1, 2, 1000);             % wide enough for Region III
xi          = omega / omega_onset;

% Continuous arrays (used for the matched-phase Nichols)
A   = zeros(size(omega));
phi = zeros(size(omega));

% Region-coloured arrays (used only for the standalone region figure).
% Pre-fill with NaN so unused indices DO NOT plot (and DO NOT collide as
% duplicate-zero frequencies inside frd()).
A1   = nan(size(omega));   phi1 = nan(size(omega));
A2   = nan(size(omega));   phi2 = nan(size(omega));
A3   = nan(size(omega));   phi3 = nan(size(omega));

m1 = xi < 1;
m2 = (xi >= 1) & (xi < 1.862);
m3 = xi >= 1.862;

for i = 1:numel(omega)
    if xi(i) < 1
        % Region I  —  no saturation
        A(i)   = 1;
        phi(i) = 0;

        A1(i)   = A(i);
        phi1(i) = phi(i);            % <-- the bug was: phi1(i) = i

    elseif xi(i) < 1.862
        % Region II —  transition (cubic in α = ω/ω_onset)
        a = xi(i);
        A(i)   = 0.2908*a^3 - 1.4396*a^2 + 1.9232*a + 0.2230;
        phi(i) = 0.5280*a^3 - 2.6213*a^2 + 3.5056*a - 1.4171;

        A2(i)   = A(i);
        phi2(i) = phi(i);

    else
        % Region III —  fully saturated  (ϖ = ω_onset/ω)
        v = omega_onset / omega(i);
        A(i)   = 4*v/pi;
        phi(i) = -acos(pi*v/2);

        A3(i)   = A(i);
        phi3(i) = phi(i);
    end
end

A_dB    = 20*log10(A);
phi_deg = rad2deg(phi);

%% =========================================================================
%  FIGURE 1 — Region-coloured Nichols-style plot of N(jω,A) alone
%             (raw axes, no phase matching — locus stays continuous)
%% =========================================================================
figure('Color','w','Name','RLE DF — region split','Position',[80 80 760 600]);
plot(rad2deg(phi1), 20*log10(A1), 'Color',[0.20 0.65 0.30], 'LineWidth', 2.2); hold on;
plot(rad2deg(phi2), 20*log10(A2), 'Color',[0.58 0.40 0.74], 'LineWidth', 2.2);
plot(rad2deg(phi3), 20*log10(A3), 'Color',[0.84 0.15 0.16], 'LineWidth', 2.2);

plot(0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
text(-3, 1.6, '\omega = \omega_{onset}', 'FontSize', 10);
phi_b = rad2deg(-acos(pi*0.537/2));
A_b   = 20*log10(4*0.537/pi);
plot(phi_b, A_b, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
text(phi_b+1, A_b+1.2, '\omega = 1.862\,\omega_{onset}', 'FontSize', 10);

grid on; box on;
xlabel('Phase  \angle N(j\omega,A)  [deg]');
ylabel('Gain  |N(j\omega,A)|  [dB]');
title('RLE describing function — three regions (standalone, no matching)');
legend({'Region I (no sat.)','Region II (transition)','Region III (full sat.)'}, ...
       'Location','southwest');

%% =========================================================================
%  FRD objects for the matched-phase Nichols
%
%  ONE continuous FRD per object — see the previous explanation.
%  KEYS:
%    • frd() takes |N|·exp(j·∠N_radians), NOT dB·exp(j·deg)
%    • frd() frequency is ω in rad/s, NOT a phase array, NOT xi
%% =========================================================================
N_frd  = frd(A .* exp(1j*phi), omega);
G_frd  = frd(Gs_ac_7, omega);
LN_frd = G_frd * N_frd;

%% =========================================================================
%  FIGURE 2 — Fig. 16 reproduction (matched-phase Nichols)
%% =========================================================================
opts = nicholsoptions;
opts.PhaseMatching       = 'on';
opts.PhaseMatchingFreq   = 1;            % rad/s
opts.PhaseMatchingValue  = -180;
opts.PhaseWrapping       = 'on';
opts.PhaseWrappingBranch = -360;

figure('Color','w','Name','Fig. 16 — q/q_c with RLE DF','Position',[80 80 900 700]);
nicholsplot(Gs_ac_7, opts);   hold on;
nicholsplot(N_frd,   opts);
nicholsplot(LN_frd,  opts);

yline(0,    '--k', 'HandleVisibility','off');
xline(-180, '--k', 'HandleVisibility','off');

% OLOP marker
plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop, 'rp', ...
     'MarkerSize', 15, 'MarkerFaceColor', 'r');
text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, ...
     sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt), ...
     'Color','r','FontWeight','bold','FontSize',12);

grid on;
xlim([-302 -62]);
ylim([-22 17]);
legend({'G_{ac\_7}', ...
        'N(j\omega,A)', ...
        'G_{ac\_7}\cdot N(j\omega,A)', ...
        'OLOP'}, 'Location','southwest');
hold off;
