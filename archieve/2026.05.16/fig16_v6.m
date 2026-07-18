clear; clc; close all;

%% =========================================================================
%  Reproduce Fig. 16 — Open-loop Nichols of q/q_c with RLE describing fn.
%  Reference: Rodrigues et al., CEAS Aeronautical Journal (2026)
%
%  v6: FIGURE 3 rewritten with manual phase matching so that G, N and G.N
%  are plotted on the SAME 360-deg branch of the Nichols chart. Each curve
%  is unwrapped independently, then shifted by an integer multiple of 360
%  so that its phase at omega_ref = 1 rad/s is closest to phi_ref = -180.
%  The OLOP marker uses the same convention as the G.N curve.
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
%% =========================================================================
omega_g = logspace(2, -1, 800);     % descending grid

A_dc_g = zeros(size(omega_g));
phi_g  = zeros(size(omega_g));
N_g    = zeros(size(omega_g));      % converged complex N(jw, A_dc)
flag_g = zeros(size(omega_g));
resn_g = zeros(size(omega_g));

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

relstep = abs(diff(A_dc_g)) ./ max(A_dc_g(1:end-1), eps);
jumps   = find(relstep > 0.5);
if ~isempty(jumps)
    fprintf('WARN: %d possible branch jumps near omega = ', numel(jumps));
    fprintf('%.3f ', omega_g(jumps));
    fprintf('rad/s\n');
end

xi_local = omega_g .* A_dc_g / R;
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
%  FIGURE 3  -  Fig. 16 reproduction with MANUAL PHASE MATCHING
%               so G, N and G.N share the same 360-deg branch.
%% =========================================================================
omega_plot = omega_g(:);

% --- Frequency responses ---
G_resp  = squeeze(freqresp(Gs_ac_7, omega_plot));   G_resp  = G_resp(:);
N_resp  = N_g(:);                                   % converged DF
GN_resp = G_resp .* N_resp;                         % G * N (open-loop)

% --- Magnitudes [dB] ---
G_mag_dB  = 20*log10(abs(G_resp));
N_mag_dB  = 20*log10(abs(N_resp));
GN_mag_dB = 20*log10(abs(GN_resp));

% --- Continuous (unwrapped) phases [deg] ---
G_phi_raw  = rad2deg(unwrap(angle(G_resp)));
N_phi_raw  = rad2deg(unwrap(angle(N_resp)));
GN_phi_raw = rad2deg(unwrap(angle(GN_resp)));

% --- Phase matching: shift each curve INDEPENDENTLY by an integer
%     multiple of 360 deg so its phase at omega_ref is nearest to phi_ref.
%     The curve shape is preserved; only the branch changes.
omega_ref = 1;        % rad/s   (common reference on the Nichols chart)
phi_ref   = -180;     % deg
[~, im]   = min(abs(omega_plot - omega_ref));

shiftN360 = @(phi_raw) phi_raw + round((phi_ref - phi_raw(im))/360)*360;

G_phi  = shiftN360(G_phi_raw);
N_phi  = shiftN360(N_phi_raw);
GN_phi = shiftN360(GN_phi_raw);

% --- OLOP marker, using the SAME convention as the G.N curve ---
olop_mag = interp1(omega_plot, GN_mag_dB, w_opt, 'pchip');
olop_phs = interp1(omega_plot, GN_phi,    w_opt, 'pchip');

fprintf('OLOP marker on G.N at omega = %.3f rad/s:\n', w_opt);
fprintf('   gain  = %+7.3f dB\n', olop_mag);
fprintf('   phase = %+7.3f deg\n\n', olop_phs);

% --- Nichols plot (manual, on top of ngrid) ---
figure('Color','w','Name','Fig. 16 - q/q_c with Gilbreath DF (phase-matched)', ...
       'Position',[80 80 950 720]);
ngrid; hold on;

hG  = plot(G_phi,  G_mag_dB,  'bo', 'MarkerSize', 4, ...
           'DisplayName', 'Aircraft Linear OLTF  G(j\omega)');
hN  = plot(N_phi,  N_mag_dB,  'r*', 'MarkerSize', 5, ...
           'DisplayName', 'Describing function  N(j\omega,A_{\delta_c})');
hGN = plot(GN_phi, GN_mag_dB, 'k.', 'MarkerSize', 8, ...
           'DisplayName', 'Aircraft Nonlinear OLTF  G\cdotN');

hOL = plot(olop_phs, olop_mag, 'rp', 'MarkerSize', 15, ...
           'MarkerFaceColor','r', ...
           'DisplayName', sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt));
text(olop_phs + 5, olop_mag, ...
     sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt), ...
     'Color','r','FontWeight','bold','FontSize',12);

yline(0,    '--', 'HandleVisibility','off');
xline(-180, '--', 'HandleVisibility','off');

xlabel('Open-Loop Phase (deg)');
ylabel('Open-Loop Gain (dB)');
title('Nichols Chart - q/q_c with Gilbreath DF (phase-matched)');
legend([hG hN hGN hOL], 'Location','best');
xlim([-300-2 -60-2]);
ylim([-20-2  15+2]);
grid on; hold off;


%% =========================================================================
%  Local functions
%% =========================================================================
function res = hb_residual(x, w, KAqc, aGc, pGc, aGa, pGa, R)
    A_dc = x(1);
    phi  = x(2);

    if A_dc <= 0
        res = [1e6; 1e6];  return;
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
    w_on = R / A_dc;
    xi   = w / w_on;
    if xi < 1
        aN = 1;            pN = 0;
    elseif xi < 1.862
        aN = 0.2908*xi^3 - 1.4396*xi^2 + 1.9232*xi + 0.2230;
        pN = 0.5280*xi^3 - 2.6213*xi^2 + 3.5056*xi - 1.4171;
    else
        v  = 1/xi;
        aN = 4*v/pi;
        pN = -acos(min(max(pi*v/2, -1), 1));
    end
end
