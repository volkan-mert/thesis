clear; clc; close all;

%% =========================================================================
%  Reproduce Fig. 16 — Open-loop Nichols of q/q_c with RLE describing fn.
%  Reference: Rodrigues et al., CEAS Aeronautical Journal (2026)
%
%  Upgrade over v3: the per-frequency describing-function amplitude A_dc
%  is now obtained from the GILBREATH iterative harmonic balance
%  (Eqs. 22–23), not from a single global omega_onset. The sweep marches
%  HIGH -> LOW omega with continuation, keeping the solver locked onto
%  the saturated branch so the Nichols "jump" is reproduced cleanly.
%% =========================================================================

%% --------- Linear plant : control law × aircraft dynamics ----------------
num_claw_7   = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7    = tf(num_claw_7, den_claw_7);

num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7  = tf(num_ldynac_7, den_ldynac_7);

Gs_ac_7      = Gs_claw_7 * Gs_ldynac_7;          % linear OLTF q / q_c

%% --------- Loop constants ------------------------------------------------
K_f  = 13.68;                 % forward (pilot) gain
A_qc = deg2rad(1);            % pilot pitch-rate command amplitude [rad/s]
R    = deg2rad(60);           % rate limit                          [rad/s]
KAqc = K_f * A_qc;            % real RHS of Eq. 22

%% --------- Linear-loop onset frequency (Duda simplified, Eq. 18) ---------
% Kept ONLY to locate the OLOP marker; the DF itself no longer uses it.
F_uc_2_urle = feedback(Gs_claw_7, Gs_ldynac_7);
fobj    = @(w) (R - KAqc*abs(evalfr(F_uc_2_urle, 1j*w))*w)^2;
fminopt = optimoptions('fmincon','Display','off');
[w_opt, ~] = fmincon(fobj, 1, [],[],[],[], 0, 100, [], fminopt);
fprintf('Closed-loop onset frequency  omega_opt = %.4f rad/s\n\n', w_opt);

%% =========================================================================
%  GILBREATH iterative harmonic balance (Eqs. 22–23)
%  --------------------------------------------------------------------
%  Two real unknowns per omega: A_dc and phi.
%     N(jw, A_dc) is the piecewise RLE DF with omega_on = R/A_dc (LOCAL).
%  Sweep direction: HIGH -> LOW.  Continuation seeds each solve from the
%  previous one, locking onto the saturated branch.
%% =========================================================================
omega_g = logspace(2, -1, 800);     % descending grid

A_dc_g = zeros(size(omega_g));
phi_g  = zeros(size(omega_g));
N_g    = zeros(size(omega_g));      % converged complex N(jw, A_dc)
flag_g = zeros(size(omega_g));
resn_g = zeros(size(omega_g));

% Initial guess at the highest frequency: Region III asymptote.
% Eq. 22 at high omega -> A_dc cos(phi) ~ K*A_qc with phi near -pi,
% so A_dc ~ K*A_qc and phi ~ -pi is a robust seed.
x0 = [KAqc; -pi];

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-14, ...
    'MaxIterations',400);

for k = 1:length(omega_g)
    w  = omega_g(k);
    Gc = evalfr(Gs_claw_7,   1j*w);
    Ga = evalfr(Gs_ldynac_7, 1j*w);

    res = @(x) hb_residual(x, w, KAqc, ...
                           abs(Gc), angle(Gc), ...
                           abs(Ga), angle(Ga), R);

    [x, fval, flag_g(k)] = fsolve(res, x0, opts);
    A_dc_g(k) = x(1);
    phi_g(k)  = x(2);
    resn_g(k) = norm(fval);

    [aN, pN]  = rle_df(w, x(1), R);
    N_g(k)    = aN * exp(1j*pN);

    % Continuation: feed this solution to the next, lower-frequency, solve
    x0 = x;
end

% Sort ascending for plotting / frd construction
[omega_g, idx] = sort(omega_g);
A_dc_g = A_dc_g(idx);
phi_g  = phi_g(idx);
N_g    = N_g(idx);
flag_g = flag_g(idx);
resn_g = resn_g(idx);

%% --------- Diagnostics ---------------------------------------------------
fprintf('Gilbreath sweep: converged @ %d / %d points (flag>0)\n', ...
        nnz(flag_g > 0), numel(flag_g));
fprintf('Max residual norm: %.2e\n', max(resn_g));

% Branch-jump detector: a large relative drop in A_dc on a descending sweep
% (= ascending after sort -> large RISE going down) typically means fsolve
% slipped to the linear branch.
relstep = abs(diff(A_dc_g)) ./ max(A_dc_g(1:end-1), eps);
jumps   = find(relstep > 0.5);
if ~isempty(jumps)
    fprintf('WARN: %d possible branch jumps near omega = ', numel(jumps));
    fprintf('%.3f ', omega_g(jumps));
    fprintf('rad/s\n');
end

% Region distribution along the converged sweep
xi_local = omega_g .* A_dc_g / R;        % xi = w / (R/A_dc) = w*A_dc/R
fprintf('Region distribution along sweep:\n');
fprintf('  I   (xi<1)        : %d points\n', nnz(xi_local <  1));
fprintf('  II  (1<=xi<1.862) : %d points\n', nnz(xi_local >= 1 & xi_local < 1.862));
fprintf('  III (xi>=1.862)   : %d points\n\n', nnz(xi_local >= 1.862));

%% =========================================================================
%  FIGURE 1  -  Bode of the converged DF
%% =========================================================================
figure('Color','w','Name','RLE DF (Gilbreath) - Bode','Position',[80 80 920 620]);

subplot(2,1,1);
semilogx(omega_g, 20*log10(abs(N_g)), 'b', 'LineWidth', 1.8); grid on; box on;
ylabel('|N(j\omega,A_{\delta_c})|  [dB]');
title('RLE describing function - Gilbreath converged solution');

subplot(2,1,2);
semilogx(omega_g, rad2deg(angle(N_g)), 'r', 'LineWidth', 1.8); grid on; box on;
xlabel('\omega  [rad/s]');
ylabel('\angle N(j\omega,A_{\delta_c})  [deg]');

%% =========================================================================
%  FIGURE 2  -  Converged unknowns A_dc(w) and phi(w)
%% =========================================================================
figure('Color','w','Name','Gilbreath converged unknowns','Position',[80 80 920 620]);

subplot(2,1,1);
semilogx(omega_g, rad2deg(A_dc_g), 'k', 'LineWidth', 1.8); grid on; box on;
ylabel('A_{\delta_c}  [deg]');
title('Converged harmonic-balance solution along sweep');

subplot(2,1,2);
semilogx(omega_g, rad2deg(phi_g), 'm', 'LineWidth', 1.8); grid on; box on;
xlabel('\omega  [rad/s]');
ylabel('\phi  [deg]');

%% =========================================================================
%  FIGURE 3  -  Fig. 16 reproduction:
%               Nichols of  G(jw),  N(jw,A),  G.N   with phase matching
%% =========================================================================
nopts = nicholsoptions;
nopts.PhaseMatching       = 'on';
nopts.PhaseMatchingFreq   = 1;             % rad/s
nopts.PhaseMatchingValue  = -180;          % anchor phase at -180 at w=1
nopts.PhaseWrapping       = 'on';
nopts.PhaseWrappingBranch = -360;

N_frd  = frd(N_g, omega_g);
G_frd  = frd(Gs_ac_7, omega_g);
LN_frd = G_frd * N_frd;                    % open-loop with Gilbreath DF

%% --- Analytic OLOP marker on the LN curve, in chart coordinates ---------
GN_vec   = squeeze(LN_frd.ResponseData(:));
phi_nat  = rad2deg(unwrap(angle(GN_vec))); % natural unwrapped phase [deg]
mag_dB   = 20*log10(abs(GN_vec));

% Apply the same phase matching the chart does: shift by integer*360
% so that phase(omega=1) = -180.
[~, im]  = min(abs(omega_g - 1));
shift    = round((-180 - phi_nat(im)) / 360) * 360;
phi_chart = phi_nat + shift;

% Then wrap to the (-360, 0] branch used by the chart options
phi_chart = mod(phi_chart, -360);

% Interpolate at w_opt
olop_mag = interp1(omega_g, mag_dB,    w_opt, 'pchip');
olop_phs = interp1(omega_g, phi_chart, w_opt, 'pchip');

fprintf('Analytic OLOP marker on G.N curve at omega = %.3f rad/s:\n', w_opt);
fprintf('   gain  = %+7.3f dB\n', olop_mag);
fprintf('   phase = %+7.3f deg\n\n', olop_phs);

%% --- Nichols plot --------------------------------------------------------
figure('Color','w','Name','Fig. 16 - q/q_c with Gilbreath DF','Position',[80 80 900 700]);
hold on;

nicholsplot(Gs_ac_7, 'bo', nopts);
legend('Aircraft Linear OLTF');

nicholsplot(N_frd,  'r*', nopts);
legend('Describing function N(j\omega,A_{\delta_c})');

nicholsplot(LN_frd, 'k.', nopts);
legend('Aircraft Nonlinear OLTF  G\cdotN');

yline(0,  '--'); legend('');
xline(-180,'--'); legend('');

plot(olop_phs, olop_mag, 'rp', 'MarkerSize', 15, 'MarkerFaceColor','r');
text(olop_phs + 5, olop_mag, ...
     sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt), ...
     'Color','r','FontWeight','bold','FontSize',12);
legend('');

grid on;
xlim([-300-2 -60-2]);
ylim([-20-2  15+2]);
hold off;


%% =========================================================================
%  Local functions
%% =========================================================================
function res = hb_residual(x, w, KAqc, aGc, pGc, aGa, pGa, R)
    % Real/imag parts of Eq. 22 / Eq. 23 of Rodrigues et al. 2026.
    A_dc = x(1);
    phi  = x(2);

    if A_dc <= 0
        res = [1e6; 1e6];  return;                 % guard the iterate
    end

    [aN, pN] = rle_df(w, A_dc, R);

    e22 = (A_dc/aGc) * cos(phi - pGc) ...
        +  A_dc * aGa * aN * cos(phi + pGa + pN) ...
        -  KAqc;
    e23 = (A_dc/aGc) * sin(phi - pGc) ...
        +  A_dc * aGa * aN * sin(phi + pGa + pN);

    res = [e22; e23];
end

function [aN, pN] = rle_df(w, A_dc, R)
    % Piecewise rate-limiter describing function (Duda / Hanke / Gilbreath).
    % Local onset frequency uses the CURRENT amplitude into the limiter.
    w_on = R / A_dc;
    xi   = w / w_on;
    if xi < 1                                      % Region I  - no sat.
        aN = 1;            pN = 0;
    elseif xi < 1.862                              % Region II - transition
        aN = 0.2908*xi^3 - 1.4396*xi^2 + 1.9232*xi + 0.2230;
        pN = 0.5280*xi^3 - 2.6213*xi^2 + 3.5056*xi - 1.4171;
    else                                           % Region III - fully sat.
        v  = 1/xi;
        aN = 4*v/pi;
        pN = -acos(min(max(pi*v/2, -1), 1));       % clamp for safety
    end
end
