clear; clc; close all;

%% =========================================================================
%  Reproduce Fig. 16 — Open-loop Nichols of q/q_c with RLE describing fn.
%  Region-coloured version: every region painted in its own colour AND
%  line style, with the locus kept continuous via a single matching offset
%  per system (G, N, G·N).
%
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
F_uc_2_urle = feedback(Gs_claw_7, Gs_ldynac_7);
R   = deg2rad(60);
Ain = 13.68/180*pi;
f       = @(w) (R - Ain*abs(evalfr(F_uc_2_urle, 1j*w))*w)^2;
fminopt = optimoptions('fmincon','Display','off');
[w_opt, ~] = fmincon(f, 5, [],[],[],[], 0, 100, [], fminopt);
fprintf('Closed-loop onset frequency  ω_onset = %.4f rad/s\n\n', w_opt);

%% =========================================================================
%  Piecewise RLE describing function   (Duda / Hanke / Gilbreath)
%% =========================================================================
omega_onset = w_opt;
omega       = logspace(-1, 2, 1000);             % must cover ω = 1 rad/s
xi          = omega / omega_onset;

A   = zeros(size(omega));
phi = zeros(size(omega));                         % radians

m1 = xi < 1;                                      % Region I
A(m1)   = 1;
phi(m1) = 0;

m2 = (xi >= 1) & (xi < 1.862);                    % Region II
a  = xi(m2);
A(m2)   = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.2230;
phi(m2) = 0.5280*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;

m3 = xi >= 1.862;                                 % Region III
v  = omega_onset ./ omega(m3);
A(m3)   = 4*v/pi;
phi(m3) = -acos(pi*v/2);

%% =========================================================================
%  Build natural complex responses on the SAME ω grid
%% =========================================================================
G_resp     = squeeze(freqresp(Gs_ac_7, omega)).';   % 1 × N
N_complex  = A .* exp(1j*phi);
LN_complex = G_resp .* N_complex;

mag_G_dB  = 20*log10(abs(G_resp));
mag_N_dB  = 20*log10(abs(N_complex));
mag_LN_dB = 20*log10(abs(LN_complex));

% Unwrapped natural phases [deg]
phi_G_nat_deg  = rad2deg(unwrap(angle(G_resp)));
phi_N_nat_deg  = rad2deg(unwrap(angle(N_complex)));
phi_LN_nat_deg = rad2deg(unwrap(angle(LN_complex)));

%% =========================================================================
%  Apply PhaseMatching @ ω = 1 rad/s to -180°,  wrap to [-360°, 0°]
%
%  KEY:  one offset per SYSTEM (G, N, LN), NOT per region.
%        That single offset is what keeps each system's locus continuous
%        across the three regions of the DF.
%% =========================================================================
match_target_deg = -180;
[~, k1] = min(abs(omega - 1));

phi_G_disp  = wrapMatch(phi_G_nat_deg,  k1, match_target_deg);
phi_N_disp  = wrapMatch(phi_N_nat_deg,  k1, match_target_deg);
phi_LN_disp = wrapMatch(phi_LN_nat_deg, k1, match_target_deg);

% Suppress vertical "snap" lines from any wrap-induced 360° jumps
phi_G_disp  = breakJumps(phi_G_disp,  180);
phi_N_disp  = breakJumps(phi_N_disp,  180);
phi_LN_disp = breakJumps(phi_LN_disp, 180);

%% =========================================================================
%  Colour & line-style per region   (edit here to taste)
%% =========================================================================
clr_R1 = [0.20 0.65 0.30];   ls_R1 = '-';     % Region I  — green,  solid
clr_R2 = [0.58 0.40 0.74];   ls_R2 = '--';    % Region II — purple, dashed
clr_R3 = [0.84 0.15 0.16];   ls_R3 = ':';     % Region III— red,    dotted
clr_G  = [0.00 0.45 0.74];                    % G_ac_7    — blue,   solid

%% =========================================================================
%  Plot
%% =========================================================================
figure('Color','w','Name','Fig. 16 — region-coloured', ...
       'Position',[80 80 940 720]);
hold on; box on; grid on;

% --- Linear plant (one continuous curve) ---
hG = plot(phi_G_disp, mag_G_dB, 'Color', clr_G, 'LineStyle','-', 'LineWidth', 2.0);

% --- N(jω,A) — three regions ---
hN1 = plot(phi_N_disp(m1), mag_N_dB(m1), 'Color',clr_R1, 'LineStyle',ls_R1, 'LineWidth',2.0);
hN2 = plot(phi_N_disp(m2), mag_N_dB(m2), 'Color',clr_R2, 'LineStyle',ls_R2, 'LineWidth',2.0);
hN3 = plot(phi_N_disp(m3), mag_N_dB(m3), 'Color',clr_R3, 'LineStyle',ls_R3, 'LineWidth',2.5);

% --- G·N(jω,A) — same hue/style per region, with markers to differentiate ---
markStride = max(1, round(numel(omega)/40));     % ~40 markers across the sweep
hL1 = plot(phi_LN_disp(m1), mag_LN_dB(m1), 'Color',clr_R1, 'LineStyle',ls_R1, ...
           'LineWidth',2.0, 'Marker','o', 'MarkerIndices', 1:markStride:nnz(m1), ...
           'MarkerFaceColor',clr_R1, 'MarkerSize',5);
hL2 = plot(phi_LN_disp(m2), mag_LN_dB(m2), 'Color',clr_R2, 'LineStyle',ls_R2, ...
           'LineWidth',2.0, 'Marker','s', 'MarkerIndices', 1:max(1,round(nnz(m2)/12)):nnz(m2), ...
           'MarkerFaceColor',clr_R2, 'MarkerSize',6);
hL3 = plot(phi_LN_disp(m3), mag_LN_dB(m3), 'Color',clr_R3, 'LineStyle',ls_R3, ...
           'LineWidth',2.5, 'Marker','d', 'MarkerIndices', 1:markStride:nnz(m3), ...
           'MarkerFaceColor',clr_R3, 'MarkerSize',6);

% --- Reference lines ---
yline(0,    '--k', 'HandleVisibility','off');
xline(-180, '--k', 'HandleVisibility','off');

% --- OLOP marker ---
hOLOP = plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop, 'rp', ...
             'MarkerSize', 16, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','k');
text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, ...
     sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt), ...
     'Color','r','FontWeight','bold','FontSize',12);

% --- Nichols grid (M-circles) ---
ngrid('on');

xlim([-302 -62]);
ylim([-22 17]);
xlabel('Open-loop phase [deg]');
ylabel('Open-loop gain [dB]');
title('Fig. 16 — q/q_c with RLE describing function (region-coloured)');

legend([hG, hN1, hN2, hN3, hL1, hL2, hL3, hOLOP], ...
       {'G_{ac\_7}', ...
        'N — Region I (no sat.)', ...
        'N — Region II (transition)', ...
        'N — Region III (full sat.)', ...
        'G_{ac\_7}\cdot N — Region I', ...
        'G_{ac\_7}\cdot N — Region II', ...
        'G_{ac\_7}\cdot N — Region III', ...
        'OLOP'}, 'Location','southwest', 'FontSize', 9);
hold off;

%% =========================================================================
%  Local helpers (require R2016b+)
%% =========================================================================
function phi_disp = wrapMatch(phi_nat_deg, k_match, target_deg)
% Phase matching at index k_match (anchor to target_deg), then wrap to [-360, 0).
    offset    = target_deg - phi_nat_deg(k_match);
    phi_match = phi_nat_deg + offset;
    phi_disp  = mod(phi_match + 360, 360) - 360;
end

function p_out = breakJumps(p, thr_deg)
% Insert NaN at indices where adjacent samples differ by more than thr_deg
% (suppresses vertical lines from wrap-induced jumps).
    p_out      = p;
    bad        = [false, abs(diff(p)) > thr_deg];
    p_out(bad) = NaN;
end
