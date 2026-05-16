%% OLOP_Nichols_Fig16.m
% Reproduces the Fig. 16 style Nichols chart for the loop
%       L(jw) = CLAW(jw) * N_rle(jw, A_i) * AC(jw)
% using Duda's piecewise rate-limiter describing function:
%   - blue   : linear open loop CLAW * AC
%   - red    : open loop with DF, CLAW * N_rle * AC
%   - green  : open-loop onset point (w = w_onset)
%
% Sweep A_i_list to also draw the OLOP locus over multiple amplitudes.

clear; clc; close all;

%% --- Plant TFs (yours) ----------------------------------------------------
num_claw_7 = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7 = [1, 21.36, 545.6, 605.7, 0];
Gs_claw   = tf(num_claw_7, den_claw_7);

num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ac     = tf(num_ldynac_7, den_ldynac_7);

L_lin = Gs_claw * Gs_ac;        % linear open loop, no RL

%% --- RL parameters --------------------------------------------------------
% Pick R and A_i so that w_onset = R/A_i lands where you want it on the
% Nichols chart (Fig. 16 in your paper used w_onset ~ 5 rad/s).
R   = 5.0;
A_i = 1.0;
w_onset = R / A_i;

% Optional: amplitudes to sweep for the OLOP locus
A_i_list = [0.5  1.0  2.0  4.0];

%% --- Frequency grid -------------------------------------------------------
% Cover well below and well above w_onset so the red curve shows the
% no-saturation -> transition -> fully saturated regions.
w = logspace(log10(w_onset/100), log10(w_onset*100), 1200);

%% --- Linear OL frequency response ----------------------------------------
H_lin = squeeze(freqresp(L_lin, w));               % complex
mag_lin_dB = 20*log10(abs(H_lin));
ph_lin_deg = unwrap(angle(H_lin)) * 180/pi;

%% --- OL with DF (red curve) ----------------------------------------------
N_rle  = duda_rl_df(w, w_onset);                   % complex DF on grid
H_DF   = H_lin .* N_rle(:);                        % L_lin * N
mag_DF_dB = 20*log10(abs(H_DF));
ph_DF_deg = unwrap(angle(H_DF)) * 180/pi;

%% --- Open-loop onset point (green) ---------------------------------------
% At w = w_onset, the cubic transition gives N ~ 1, so OLOP sits on the
% linear curve. Evaluate exactly there to avoid grid-snapping.
H_onset_lin = squeeze(freqresp(L_lin, w_onset));
N_onset     = duda_rl_df(w_onset, w_onset);
H_onset_DF  = H_onset_lin * N_onset;

mag_OLOP_dB = 20*log10(abs(H_onset_DF));
ph_OLOP_deg = angle(H_onset_DF) * 180/pi;
% Match unwrapped branch of the linear curve so the marker lands on it
ph_OLOP_deg = wrap_to_branch(ph_OLOP_deg, ph_lin_deg, w, w_onset);

%% --- OLOP locus over a range of input amplitudes -------------------------
N_sweep   = numel(A_i_list);
OLOP_pts  = zeros(N_sweep, 2);                     % [phase_deg, mag_dB]
for k = 1:N_sweep
    Ai_k     = A_i_list(k);
    won_k    = R / Ai_k;
    Hk_lin   = squeeze(freqresp(L_lin, won_k));
    Nk_onset = duda_rl_df(won_k, won_k);
    Hk       = Hk_lin * Nk_onset;
    ph_k     = wrap_to_branch(angle(Hk)*180/pi, ph_lin_deg, w, won_k);
    OLOP_pts(k,:) = [ph_k, 20*log10(abs(Hk))];
end

%% --- Plot -----------------------------------------------------------------
figure('Position',[100 100 900 650]); hold on; grid on; box on;

plot(ph_lin_deg, mag_lin_dB, '-o', ...
    'Color',[0 0.45 0.85], 'MarkerSize',4, 'LineWidth',1.0, ...
    'MarkerFaceColor','none', 'DisplayName','Linear frequency response');

plot(ph_DF_deg, mag_DF_dB, '-.', ...
    'Color',[0.85 0.1 0.1], 'LineWidth',1.6, ...
    'DisplayName','Describing function');
plot(ph_DF_deg(1:15:end), mag_DF_dB(1:15:end), 'o', ...
    'Color',[0.85 0.1 0.1], 'MarkerSize',4, 'MarkerFaceColor',[0.85 0.1 0.1], ...
    'HandleVisibility','off');

plot(ph_OLOP_deg, mag_OLOP_dB, 'o', ...
    'MarkerSize',10, 'MarkerFaceColor',[0.1 0.75 0.2], ...
    'MarkerEdgeColor','k', 'LineWidth',1.0, ...
    'DisplayName','Open-loop onset point');

if numel(A_i_list) > 1
    plot(OLOP_pts(:,1), OLOP_pts(:,2), 's--', ...
        'Color',[0.1 0.5 0.1], 'MarkerSize',7, ...
        'MarkerFaceColor',[0.6 0.95 0.6], ...
        'DisplayName','OLOP locus (A_i sweep)');
end

% Reference lines
yl = ylim; xl = xlim;
plot([-180 -180], yl, 'k:', 'HandleVisibility','off');
plot(xl, [0 0],   'k:',   'HandleVisibility','off');

xlabel('Open-loop phase (\circ)');
ylabel('Open-loop gain (dB)');
title(sprintf(['Nichols chart: q/q_c   |   ' ...
               'A_i = %.2f,  R = %.2f,  \\omega_{onset} = %.3f rad/s'], ...
              A_i, R, w_onset));
legend('Location','southwest');

% Optional Nichols M-circles. Uncomment if you want them:
% hN = nichols(L_lin, w); ngrid;  % spawns a separate figure with the grid

%% --- Local functions -----------------------------------------------------
function [N, A, phi] = duda_rl_df(w, w_onset)
% Piecewise describing function of an ideal rate limiter (Duda).
%   alpha = w/w_onset = 1/wbar
    alpha = w(:).' / w_onset;
    A   = zeros(size(alpha));
    phi = zeros(size(alpha));

    idx1 = alpha < 1;                                    % no saturation
    A(idx1)   = 1;
    phi(idx1) = 0;

    idx2 = (alpha >= 1) & (alpha < 1.862);               % cubic transition
    a = alpha(idx2);
    A(idx2)   = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.2230;
    phi(idx2) = 0.5280*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;

    idx3 = alpha >= 1.862;                               % full saturation
    wbar = 1 ./ alpha(idx3);
    A(idx3)   = 4*wbar/pi;
    phi(idx3) = -acos(pi*wbar/2);

    N = A .* exp(1j*phi);
end

function ph_branch = wrap_to_branch(ph_deg, ph_ref_unwrapped, w_ref, w_eval)
% Place a single phase value on the same 360 deg branch as the unwrapped
% reference curve at frequency w_eval, by snapping to the nearest sample.
    [~,i0] = min(abs(w_ref - w_eval));
    ph_branch = ph_deg + 360*round((ph_ref_unwrapped(i0) - ph_deg)/360);
end