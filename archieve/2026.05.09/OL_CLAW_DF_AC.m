%% OL_CLAW_DF_AC.m
% Open-loop CLAW * N_rle * AC using Duda's piecewise rate-limiter
% describing function (DF). The DF is amplitude-dependent, so it is
% built as an FRD on the same frequency grid as CLAW and AC.
%
% Regions (alpha = w/w_onset, w_onset = R/A_i):
%   I  : alpha < 1                    -> A=1, phi=0          (no saturation)
%   II : 1 <= alpha < 1.862           -> cubic spline         (transition)
%   III: alpha >= 1.862               -> A=4*wbar/pi,
%                                        phi=-acos(pi*wbar/2)  (fully saturated)
% with wbar = w_onset/w = 1/alpha.
%
% Reference: Duda, "Effects of Rate Limiting Elements in Flight Control
% Systems - A New PIO-Criterion", AIAA-95-3304; AGARD-AR-335.

clear; clc; close all;

%% --- Plant transfer functions (your CLAW and AC) ---
num_claw_7 = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7 = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7  = tf(num_claw_7, den_claw_7);

num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7  = tf(num_ldynac_7, den_ldynac_7);

%% --- Rate-limiter parameters ---
R   = 1.0;     % rate limit [units of u_rle per second]
A_i = 1.0;     % sinusoidal input amplitude into the RL  -> sweep for OLOP
w_onset = R / A_i;

%% --- Frequency grid ---
w = logspace(-1, 2, 1000);   % rad/s

%% --- Piecewise describing function N_rle(jw, A_i) ---
[N_rle, A_mag, phi_rad] = duda_rl_df(w, w_onset);
DF_frd = frd(N_rle, w);

%% --- FRDs of CLAW and AC on the same grid, then product ---
CLAW_frd = frd(Gs_claw_7,   w);
AC_frd   = frd(Gs_ldynac_7, w);

OL_no_RL = CLAW_frd * AC_frd;            % linear OL (no rate limiter)
OL_with  = CLAW_frd * DF_frd * AC_frd;   % OL with describing function

%%
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

%% --- Plots ---
figure('Name','Bode'); 
bode(OL_no_RL, 'b--', OL_with, 'r-', w); grid on;
legend('CLAW \cdot AC','CLAW \cdot N_{rle} \cdot AC','Location','best');
title(sprintf('Bode | A_i=%.3g, R=%.3g, \\omega_{onset}=%.3g rad/s', ...
              A_i, R, w_onset));
legend show
xlim auto
ylim auto

figure('Name','Nichols');
opts = nicholsoptions;
opts.PhaseMatching      = 'on';
opts.PhaseMatchingFreq  = 1;        % rad/s
opts.PhaseMatchingValue = -180;     % anchor phase at that frequency
opts.PhaseWrapping = 'on';
opts.PhaseWrappingBranch = -360; % Wraps between -360 and 0
nicholsplot(OL_no_RL, 'b--', opts);
hold on;
nicholsplot(OL_with, 'r-', opts);
% nichols(OL_no_RL, 'b--', OL_with, 'r-', w);
%
plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop,'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
hold on
text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, sprintf('\\omega_{OLOP} = %.3f', w_opt), 'Color', 'red', 'FontWeight', 'bold', 'FontSize', 12);
%
grid on; 
ngrid;
legend('CLAW \cdot AC','CLAW \cdot N_{rle} \cdot AC','Location','best');
title('Nichols');
legend show
%
xlim([-300 -50]);
ylim([-20 15]);
%


figure('Name','Nyquist');
nyquist(OL_no_RL, 'b--', OL_with, 'r-', w); grid on;
legend('CLAW \cdot AC','CLAW \cdot N_{rle} \cdot AC','Location','best');
title('Nyquist');
legend show
xlim auto
ylim auto


%% --- OLOP markers (open-loop onset point) ---
% OLOP convention: plot |CLAW*AC|, angle(CLAW*AC * N_rle) at the
% closed-loop neutral-stability frequency / onset frequency.
% Here we just mark the response at w_onset and at the linear -180 deg
% crossing for quick inspection.
[mag180, ph180, w180] = margin(Gs_claw_7 * Gs_ldynac_7);  

if isfinite(w180)
    H_at_w180  = squeeze(freqresp(OL_with,  w180));
    H_lin_w180 = squeeze(freqresp(OL_no_RL, w180));
    fprintf('At linear w_180 = %.4f rad/s:\n', w180);
    fprintf('  |CLAW*AC|       = %7.3f dB, ang = %7.2f deg\n', ...
            20*log10(abs(H_lin_w180)), rad2deg(angle(H_lin_w180)));
    fprintf('  |CLAW*N*AC|     = %7.3f dB, ang = %7.2f deg\n', ...
            20*log10(abs(H_at_w180)), rad2deg(angle(H_at_w180)));
end

H_at_wons = squeeze(freqresp(OL_with, w_onset));
fprintf('At w_onset = %.4f rad/s:\n', w_onset);
fprintf('  |CLAW*N*AC|     = %7.3f dB, ang = %7.2f deg\n', ...
        20*log10(abs(H_at_wons)), rad2deg(angle(H_at_wons)));

%% ===== Local function =====
function [N, A, phi] = duda_rl_df(w, w_onset)
% DUDA_RL_DF  Piecewise describing function of an ideal rate limiter.
%   N   : complex DF, N = A .* exp(1j*phi)
%   A   : magnitude
%   phi : phase in radians
%   w     : frequency vector (rad/s)
%   w_onset = R / A_i
%
% Cubic-spline transition coefficients per Duda (Region II).

    alpha = w(:).' ./ w_onset;     % alpha = w / w_onset (= 1/wbar)

    A   = zeros(size(alpha));
    phi = zeros(size(alpha));

    % --- Region I: no saturation
    idx1 = alpha < 1;
    A(idx1)   = 1;
    phi(idx1) = 0;

    % --- Region II: cubic spline in alpha
    idx2 = (alpha >= 1) & (alpha < 1.862);
    a    = alpha(idx2);
    A(idx2)   = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.2230;
    phi(idx2) = 0.5280*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;   % rad

    % --- Region III: fully developed saturation
    idx3 = alpha >= 1.862;
    wbar = 1 ./ alpha(idx3);                 % wbar = w_onset/w
    A(idx3)   = 4*wbar / pi;
    phi(idx3) = -acos(pi*wbar/2);

    N = A .* exp(1j*phi);
end
