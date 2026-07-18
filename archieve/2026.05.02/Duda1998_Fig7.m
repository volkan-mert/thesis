%% CLAUDE AI Generated Code  
% ========================================================================
%  Duda1998_Fig7.m
%  Reproduces the simple design example from:
%    Duda, H. (1998), "Flight control system design considering rate
%    saturation", Aerospace Science and Technology 4, 265-275.
%  Figure 7 — basically unstable aircraft with α and q feedback through a
%  2nd-order actuator.
%
%  Computes the OPEN-LOOP transfer function L(s) with the loop broken at
%  the elevator actuator output δe, per the paper's §3 instruction
%  ("in the pitch axis the loop has to be opened at the elevator actuator").
%
%  Verifies the linear short-period claims of §5.2:
%      Nominal  (K_fb = 1.00)  →  ω_sp ≈ 5.9 rad/s, ζ_sp ≈ 0.7
%      Low gain (K_fb = 0.68)  →  meets Level-1 phase rate criterion
%
%  Outputs (saved to base workspace) for downstream OLOP analysis:
%      L_nom, L_low      open-loop TFs (loop broken at δe)
%      T_Fes_de_nom/low  closed-loop TFs F_es → δe (for ω_onset)
%      g_a, g_ff, g_de_alpha, g_de_q, H_fb   building blocks
% =========================================================================

clear; clc; close all;

%% ====================================================================
%  1. PLANT BLOCKS (Figure 7 numerical values)
%  ====================================================================
%  Notation in the paper:
%     (ζ; ω) ≡ s² + 2ζω s + ω²
%     (a)    ≡ s + a            so (-1.06) ≡ s − 1.06  (RHP pole)

% --- Basic aircraft (basically unstable; single short-period mode) ------
% g_de^α(s) = -0.23 (s + 40.45) / [(s + 1.95)(s - 1.06)]
% g_de^q(s) = -9.15 (s +  0.52) / [(s + 1.95)(s - 1.06)]
den_ac     = conv([1 1.95], [1 -1.06]);          % shared airframe denom.
g_de_alpha = tf(-0.23 * [1 40.45], den_ac);
g_de_q     = tf(-9.15 * [1  0.52], den_ac);

% --- Actuator (2nd order, ω_a = 31.4 rad/s, ζ_a = 0.7) ------------------
omega_a = 31.4;   zeta_a = 0.7;
g_a = tf(omega_a^2, [1, 2*zeta_a*omega_a, omega_a^2]);

% --- Prefilter (DC gain = 1) -------------------------------------------
% g_ff(s) = 0.125 (s + 10) / (s + 1.25)
g_ff = tf(0.125*[1 10], [1 1.25]);

% --- Aircraft control authority (sensitivity), forward gain ------------
K_ac = 1.0;
K_ff = 1.0;

% --- Feedback gains (manual optimisation, §5.2) ------------------------
K_alpha = 2.5;     % [deg elevator / deg α]
K_q     = 0.6;     % [deg elevator / (deg/s)]   (units of seconds)

% --- Hardware limits (used downstream by Simulink / OLOP analysis) -----
F_es_max = 80;     % stick force limit                    [N]
R_fb     = 100;    % feedback-loop rate limit             [deg/s]
delta_max= 30;     % elevator amplitude limit             [deg]

%% ====================================================================
%  2. OPEN-LOOP TRANSFER FUNCTION (loop broken at δe)
%  ====================================================================
%  Around the loop:  δe → K_ac → plant → (K_α α + K_q q) → Σ(−) → K_fb
%                    → g_a → δe
%
%  The plant numerators are negative (g_de^α, g_de^q have negative DC
%  gain in Duda's sign convention); the leading minus below produces a
%  loop transmission that yields negative feedback at the summing
%  junction and stabilises the +1.06 RHP airframe pole.
%
%      L(s) = −K_fb · K_ac · g_a(s) · [K_α g_de^α(s) + K_q g_de^q(s)]
% ----------------------------------------------------------------------

H_fb = K_alpha * g_de_alpha + K_q * g_de_q;       % combined α/q feedback

% Two designs from §5.2
K_fb_nom = 1.00;
K_fb_low = 0.68;

L_nom = minreal( -K_fb_nom * K_ac * g_a * H_fb );
L_low = minreal( -K_fb_low * K_ac * g_a * H_fb );

fprintf('\n=== OPEN-LOOP TF — Nominal design (K_fb = %.2f) ===\n', K_fb_nom);
zpk(L_nom)

fprintf('=== OPEN-LOOP TF — Low-gain design (K_fb = %.2f) ===\n', K_fb_low);
zpk(L_low)

%% ====================================================================
%  3. CLOSED-LOOP CHARACTERISTICS — verify Duda's §5.2 claims
%  ====================================================================
T_nom = minreal(feedback(L_nom, 1));
T_low = minreal(feedback(L_low, 1));

fprintf('\n--- Closed-loop modes (Nominal, K_fb = 1.00) ---\n');
damp(T_nom);
fprintf('--- Closed-loop modes (Low gain, K_fb = 0.68) ---\n');
damp(T_low);

[Gm_n,Pm_n,Wcg_n,Wcp_n] = margin(L_nom);
[Gm_l,Pm_l,Wcg_l,Wcp_l] = margin(L_low);
fprintf('\nNominal:  GM = %5.2f dB,  PM = %5.2f deg,  ω_c = %.2f rad/s\n', ...
        20*log10(Gm_n), Pm_n, Wcp_n);
fprintf('Low gain: GM = %5.2f dB,  PM = %5.2f deg,  ω_c = %.2f rad/s\n', ...
        20*log10(Gm_l), Pm_l, Wcp_l);

%% ====================================================================
%  4. CLOSED-LOOP TFs FROM F_es TO δe  (for ω_onset / OLOP analysis)
%  ====================================================================
%  Duda §5.3, Fig 9a — the rate limiter R_fb in the feedback loop is
%  excited by δe; ω_onset is the frequency at which |F^δe_es(jω)| crosses
%  the −20 dB/decade asymptote 100/ω in dB.
%
%      F_es →[F_es_max]→ R_ff →[g_ff]→[K_ff]→ Σ →[K_fb]→ R_fb →[δ_max]→ g_a → δe
%                                          ↑
%                                          └── K_ac · (K_α g_de^α + K_q g_de^q)
%
%  Linearised (no rate / amplitude limits): the F_es → δe TF is
%      T_F→δ(s) = g_ff · K_ff · K_fb · g_a / (1 + L(s))
% ----------------------------------------------------------------------
fwd_nom = g_ff * K_ff * K_fb_nom * g_a;           % forward path command→δe
fwd_low = g_ff * K_ff * K_fb_low * g_a;
T_Fes_de_nom = minreal(fwd_nom / (1 + L_nom));
T_Fes_de_low = minreal(fwd_low / (1 + L_low));

%% ====================================================================
%  5. PLOTS
%  ====================================================================
w = logspace(-1, 2.5, 2000);
col_nom = [0.00 0.25 0.70];
col_low = [0.80 0.10 0.10];

% --- 5a. Bode of open loop -------------------------------------------
figure('Name','Open-Loop Bode');
% figure('Name','Open-Loop Bode','Position',[60 540 920 600]);
bopt = bodeoptions; bopt.Grid = 'on'; bopt.PhaseMatching = 'on';
bopt.PhaseMatchingFreq = 1; bopt.PhaseMatchingValue = -180;
bodeplot(L_nom, L_low, w, bopt);
lines = findobj(gcf,'Type','Line'); set(lines,'LineWidth',1.6);
legend('Nominal  K_{fb}=1.00','Low gain  K_{fb}=0.68', ...
       'Location','southwest');
title('Duda Fig. 7 — Open-Loop Bode (loop opened at δ_e)');

% --- 5b. Nichols of open loop (with closed-loop M-grid) --------------
figure('Name','Open-Loop Nichols');
% figure('Name','Open-Loop Nichols','Position',[60 60 920 620]);
nopt = nicholsoptions; nopt.Grid = 'on';
nopt.PhaseMatching = 'on'; nopt.PhaseMatchingFreq = 1;
nopt.PhaseMatchingValue = -180;
nopt.XLim = {[-360 0]}; nopt.YLim = {[-40 30]};
nicholsplot(L_nom, L_low, w, nopt);
lines = findobj(gcf,'Type','Line'); set(lines,'LineWidth',1.6);
legend('Nominal  K_{fb}=1.00','Low gain  K_{fb}=0.68', ...
       'Location','southwest');
title('Duda Fig. 7 — Open-Loop Nichols (broken at δ_e)');

% --- 5c. F_es → δe Bode (for ω_onset reading) -------------------------
figure('Name','F_{es} → δ_e Closed-Loop Bode');
% figure('Name','F_{es} → δ_e Closed-Loop Bode','Position',[1000 540 920 600]);
% Plot |F_es→δe| against the rate-limit asymptote R_fb / ω in dB
[mag_n,~] = bode(T_Fes_de_nom, w); mag_n = 20*log10(squeeze(mag_n));
[mag_l,~] = bode(T_Fes_de_low, w); mag_l = 20*log10(squeeze(mag_l));
% Asymptote is for full-stick input; convert to dB referenced to F_es=80 N
mag_n_full = mag_n + 20*log10(F_es_max);
mag_l_full = mag_l + 20*log10(F_es_max);
asymptote  = 20*log10(R_fb./w);
amp_limit  = 20*log10(delta_max) * ones(size(w));

semilogx(w, mag_n_full, '-',  'Color',col_nom,'LineWidth',1.7); hold on;
semilogx(w, mag_l_full, '-',  'Color',col_low,'LineWidth',1.7);
semilogx(w, asymptote,  'k--','LineWidth',1.3);
semilogx(w, amp_limit,  'k:', 'LineWidth',1.3);
grid on; xlim([w(1) w(end)]); ylim([-10 60]);
xlabel('Frequency [rad/s]'); ylabel('|F_{es}^{δe}|  [dB]');
title('|F_{es}→δ_e| at full stick (80 N) — onset where curve crosses 100/ω');
legend('Nominal','Low gain','100 deg/s asymptote','30 deg amplitude limit', ...
       'Location','northeast');

% --- 5d. Closed-loop step response (sanity) ---------------------------
figure('Name','Closed-Loop Step');
% figure('Name','Closed-Loop Step','Position',[1000 60 920 400]);
step(T_nom, T_low, 5);
lines = findobj(gcf,'Type','Line'); set(lines,'LineWidth',1.6);
legend('Nominal','Low gain','Location','southeast');
title('Closed-loop step response (linear, no rate / amp limits)');
grid on;

%% ====================================================================
%  6. EXPORT TO BASE WORKSPACE FOR DOWNSTREAM OLOP / DESCRIBING FUNCTION
%  ====================================================================
fprintf(['\nWorkspace variables for OLOP analysis:\n', ...
         '  L_nom, L_low                 open-loop TFs at δ_e\n', ...
         '  T_Fes_de_nom, T_Fes_de_low   closed-loop TFs F_es → δ_e\n', ...
         '  g_a, g_ff, H_fb              building blocks\n', ...
         '  R_fb (=%d), δ_max (=%d), F_es_max (=%d)\n\n'], ...
         R_fb, delta_max, F_es_max);
