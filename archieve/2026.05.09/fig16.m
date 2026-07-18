clear; close all; clc;

%% =========================================================================
%  Reproduce Fig. 16 — Open-loop Nichols of q/q_c with RLE describing fn.
%  Reference: Rodrigues et al., CEAS Aeronautical Journal (2026).
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
F_uc_2_urle = feedback(Gs_claw_7, Gs_ldynac_7);  % from u_c to u_RLE

R   = deg2rad(60);            % rate limit               [rad/s]
Ain = 13.68/180*pi;           % pilot input × K_f        [rad]

f       = @(w) (R - Ain*abs(evalfr(F_uc_2_urle, 1j*w))*w)^2;
fminopt = optimoptions('fmincon','Display','off');
[w_opt, ~] = fmincon(f, 5, [],[],[],[], 0, 100, [], fminopt);

fprintf('Closed-loop onset frequency  ω_onset = %.4f rad/s\n\n', w_opt);

%% =========================================================================
%  Piecewise RLE describing function   (vectorized — no sparse-array trap)
%% =========================================================================
omega_onset = w_opt;                             % from Eq. 18
omega       = logspace(-1, 2, 1000);             % must cover ω = 1 rad/s
xi          = omega / omega_onset;               % ω / ω_onset

A   = zeros(size(omega));                        % |N(jω,A)|   linear
phi = zeros(size(omega));                        % ∠N(jω,A)    radians

% Region I  —  no saturation
m1      = xi < 1;
A(m1)   = 1;
phi(m1) = 0;

% Region II —  transition (cubic in α = ω/ω_onset)
m2 = (xi >= 1) & (xi < 1.862);
a  = xi(m2);
A(m2)   = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.2230;
phi(m2) = 0.5280*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;

% Region III —  fully saturated  (ϖ = ω_onset/ω)
m3 = xi >= 1.862;
v  = omega_onset ./ omega(m3);
A(m3)   = 4*v/pi;
phi(m3) = -acos(pi*v/2);

A_dB    = 20*log10(A);
phi_deg = rad2deg(phi);



%% =========================================================================
%  FIGURE 1 — Standalone Nichols-style locus (region-coloured)
%             Use THIS figure to see the three regions of N(jω,A).
%             Phase matching is *not* applied here, so the locus is
%             continuous across regions.
%% =========================================================================
figure('Color','w','Name','RLE DF — Nichols-style','Position',[80 80 760 600]);

plot(phi_deg(m1), A_dB(m1), 'Color',[0.20 0.65 0.30], 'LineWidth', 2.2); hold on;
plot(phi_deg(m2), A_dB(m2), 'Color',[0.58 0.40 0.74], 'LineWidth', 2.2); hold on;
plot(phi_deg(m3), A_dB(m3), 'Color',[0.84 0.15 0.16], 'LineWidth', 2.2); hold on;

plot(0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6); hold on;
text(-3, 1.6, '\omega = \omega_{onset}', 'FontSize', 10); hold on;

phi_b = rad2deg(-acos(pi*0.537/2));
A_b   = 20*log10(4*0.537/pi);
plot(phi_b, A_b, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6); hold on;
text(phi_b+1, A_b+1.2, '\omega = 1.862\,\omega_{onset}', 'FontSize', 10); hold on;

grid on; box on;
xlabel('Phase  \angle N(j\omega,A)  [deg]');
ylabel('Gain  |N(j\omega,A)|  [dB]');
title('RLE describing function — three regions on Nichols-style axes');
legend({'Region I (no sat.)', ...
        'Region II (transition)', ...
        'Region III (full sat.)'}, 'Location','southwest');
legend show

%% =========================================================================
%  FRD objects for the matched-phase Nichols
%
%  Why ONE FRD per object instead of three (per region):
%    PhaseMatching anchors phase at ω = 1 rad/s. That frequency only lies
%    inside Region I's range. For per-region FRDs of Region II / III,
%    MATLAB falls back to the closest in-range frequency, which gives
%    each region a DIFFERENT matching offset — the locus then breaks
%    into three disjoint segments on the matched chart. A single FRD
%    over the full ω grid keeps the locus continuous.
%
%  Two more landmines that bit the previous draft:
%    (a) frd() expects |N| · exp(j·∠N_radians) — NOT  |N|_dB · exp(j·deg)
%    (b) frd() frequency must be ω in rad/s, NOT the dimensionless
%        xi = ω/ω_onset (otherwise multiplication with G_frd is junk).
%% =========================================================================
N_frd  = frd(A .* exp(1j*phi), omega);           % describing function
G_frd  = frd(Gs_ac_7, omega);                    % plant on same grid
LN_frd = G_frd * N_frd;                          % open-loop with RLE DF

%% =========================================================================
%  FIGURE 2 — Fig. 16 reproduction (matched-phase Nichols)
%% =========================================================================
opts = nicholsoptions;
opts.PhaseMatching       = 'on';
opts.PhaseMatchingFreq   = 1;            % rad/s
opts.PhaseMatchingValue  = -180;         % anchor at -180° at ω = 1 rad/s
opts.PhaseWrapping       = 'on';
opts.PhaseWrappingBranch = -360;

figure('Color','w','Name','Fig. 16 — q/q_c with RLE DF','Position',[80 80 900 700]);
nicholsplot(Gs_ac_7, opts);   hold on;
nicholsplot(N_frd,   opts);   hold on;
nicholsplot(LN_frd,  opts);   hold on;

yline(0,    '--k', 'HandleVisibility','off');   hold on;
xline(-180, '--k', 'HandleVisibility','off');   hold on;

% OLOP marker (cursor-picked on the published, matched-phase chart)
plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop, 'rp', ...
     'MarkerSize', 15, 'MarkerFaceColor', 'r');    hold on;
text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, ...
     sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt), ...
     'Color','r','FontWeight','bold','FontSize',12);   hold on;

grid on;
xlim([-302 -62]);
ylim([-22 17]);
legend({'G_{ac\_7}', ...
        'N(j\omega,A)', ...
        'G_{ac\_7}\cdot N(j\omega,A)', ...
        'OLOP'}, 'Location','southwest');
hold off;

%% =========================================================================
%  FIGURE 3 — Region-split Open-Loop Nichols
%  (Replica of Figure 2, but LN locus is colored by RLE region)
%% =========================================================================
figure('Color','w','Name','Figure 3 — Region split','Position',[80 80 900 700]);

% Draw the plant and DF alone as background reference
nicholsplot(G_frd, opts);   
hold on;
nicholsplot(N_frd,   opts);

% Manually compute matched phase for LN_frd to avoid nicholsplot breaking 
% the locus when the FRD is partitioned.
LN_resp = squeeze(LN_frd.ResponseData).';
LN_mag_dB = 20*log10(abs(LN_resp));
LN_phase_raw = rad2deg(unwrap(angle(LN_resp)));

% Phase matching at w = 1 rad/s
[~, idx1] = min(abs(omega - 1));
phase_shift = -180 - LN_phase_raw(idx1);
LN_phase_matched = LN_phase_raw + phase_shift;

% Wrap phase specifically to [-360, 0) branch
LN_phase_wrapped = mod(LN_phase_matched, 360) - 360;

% Plot the 3 regions with custom markers. 
% Added MarkerIndices to prevent marker clutter (1000 points is too many for every point)
plot(LN_phase_wrapped(m1), LN_mag_dB(m1), 'r-->', 'LineWidth', 2.2, 'MarkerIndices', 1:15:sum(m1));
plot(LN_phase_wrapped(m2), LN_mag_dB(m2), 'g--*', 'LineWidth', 2.2, 'MarkerIndices', 1:5:sum(m2));
plot(LN_phase_wrapped(m3), LN_mag_dB(m3), 'b--o', 'LineWidth', 2.2, 'MarkerIndices', 1:15:sum(m3));

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

% Create exact dummy lines for a flawless legend that avoids the Warning
h_leg(1) = plot(NaN, NaN, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); % Default blue for G
h_leg(2) = plot(NaN, NaN, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5); % Default orange for N
h_leg(3) = plot(NaN, NaN, 'r-->', 'LineWidth', 2.2);
h_leg(4) = plot(NaN, NaN, 'g--*', 'LineWidth', 2.2);
h_leg(5) = plot(NaN, NaN, 'b--o', 'LineWidth', 2.2);
h_leg(6) = plot(NaN, NaN, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

legend(h_leg, {'G_{ac\_7}', ...
               'N(j\omega,A)', ...
               'G_{ac\_7}\cdot N (Region I)', ...
               'G_{ac\_7}\cdot N (Region II)', ...
               'G_{ac\_7}\cdot N (Region III)', ...
               'OLOP'}, 'Location','southwest');
hold off;


%% =========================================================================
%  FIGURE 4 — Fig. 16 reproduction (matched-phase Nichols) - SPLIT
%% =========================================================================
opts = nicholsoptions;
opts.PhaseMatching       = 'on';
opts.PhaseMatchingFreq   = 1;            % rad/s
opts.PhaseMatchingValue  = -180;         % anchor at -180° at ω = 1 rad/s
opts.PhaseWrapping       = 'on';
opts.PhaseWrappingBranch = -360;

figure('Color','w','Name','Fig. 16 — q/q_c with RLE DF','Position',[80 80 900 700]);
nicholsplot(G_frd, opts);   hold on;
nicholsplot(N_frd,   opts);   hold on;
nicholsplot(LN_frd,  opts);   hold on;

yline(0,    '--k', 'HandleVisibility','off');   hold on;
xline(-180, '--k', 'HandleVisibility','off');   hold on;

% OLOP marker (cursor-picked on the published, matched-phase chart)
plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop, 'rp', ...
     'MarkerSize', 15, 'MarkerFaceColor', 'r');    hold on;
text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, ...
     sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt), ...
     'Color','r','FontWeight','bold','FontSize',12);   hold on;

grid on;
xlim([-302 -62]);
ylim([-22 17]);
legend({'G_{ac\_7}', ...
        'N(j\omega,A)', ...
        'G_{ac\_7}\cdot N(j\omega,A)', ...
        'OLOP'}, 'Location','southwest');
hold off;
