clear; clc; close all;

%% =========================================================================
%  Reproduce Fig. 16 — Open-loop Nichols of q/q_c with RLE describing fn.
%  Reference: Rodrigues et al., CEAS Aeronautical Journal (2026),
%             "Conditional integrator sliding mode control to reduce
%              susceptibility to pilot-induced oscillations"
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
Ain = 13.68/180*pi;           % pilot input amplitude × K_f   [rad]

f       = @(w) (R - Ain*abs(evalfr(F_uc_2_urle, 1j*w))*w)^2;
fminopt = optimoptions('fmincon','Display','off');
[w_opt, ~] = fmincon(f, 1, [],[],[],[], 0, 100, [], fminopt);

fprintf('Closed-loop onset frequency  ω_onset = %.4f rad/s\n\n', w_opt);

%% =========================================================================
%  Harmonic Balance Numerical Solver (Eqs 22 and 23) via Anonymous Functions
%% =========================================================================
omega_onset = w_opt;                             
omega       = logspace(-1, 2, 1000);              

% 1. Map parameters from the plant to the balance equations
K_gain = 13.68;              % Forward gain K_f
A_qc   = 1.0 * (pi/180);     % Pilot pitch rate command amplitude (1 deg -> rad)
R_rate = 60.0 * (pi/180);    % Maximum rate limit (deg/s -> rad/s)

% Preallocate arrays
A_delta_c = zeros(size(omega));
phi_sol   = zeros(size(omega));
A         = zeros(size(omega)); % |N(jω,A)| of the RLE
phi       = zeros(size(omega)); % ∠N(jω,A) of the RLE

% 2. Anonymous function for the piecewise Describing Function N(Adc)
% Evaluates Region I, II, and III using boolean masking. Returns [mag_N; phi_N]
calc_N = @(Adc, w, RateL, ratio) ...
    (ratio < 1) .* [1; 0] + ...
    (ratio >= 1 & ratio < 1.862) .* ...
        [0.2908*ratio^3 - 1.4396*ratio^2 + 1.9232*ratio + 0.2230; ...
         0.5280*ratio^3 - 2.6213*ratio^2 + 3.5056*ratio - 1.4171] + ...
    (ratio >= 1.862) .* ...
        [4*(1./ratio)/pi; -acos(min(1, pi*(1./ratio)/2))]; 

% 3. Anonymous function for the objective (Eqs 22 & 23)
% Uses feval to evaluate calc_N once and distribute its output [mag_N; phi_N]
harmonic_eqs = @(X, w, mGc, pGc, mGa, pGa) feval(...
    @(Adc, ph, N_vals) [ ...
        (Adc / mGc) * cos(ph - pGc) + Adc * mGa * N_vals(1) * cos(ph + pGa + N_vals(2)) - K_gain * A_qc; ...
        (Adc / mGc) * sin(ph - pGc) + Adc * mGa * N_vals(1) * sin(ph + pGa + N_vals(2)) ...
    ], X(1), X(2), calc_N(X(1), w, R_rate, w * X(1) / R_rate));

% 4. Solve the coupled non-linear equations
% Sweeping downwards to catch the saturated branch naturally
X0 = [R_rate / omega(end), -pi/2]; 
opts_fsolve = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'trust-region-dogleg');

for i = length(omega):-1:1
    w = omega(i);
    
    % Evaluate linear plant and controller at current frequency
    Gc_val = evalfr(Gs_claw_7, 1j*w);
    Ga_val = evalfr(Gs_ldynac_7, 1j*w);
    
    mGc = abs(Gc_val); pGc = angle(Gc_val);
    mGa = abs(Ga_val); pGa = angle(Ga_val);
    
    % Solve for X = [A_delta_c; phi_sol]
    [X_sol, ~, exitflag] = fsolve(@(X) harmonic_eqs(X, w, mGc, pGc, mGa, pGa), X0, opts_fsolve);
    
    if exitflag > 0
        A_delta_c(i) = X_sol(1);
        phi_sol(i)   = X_sol(2);
        X0 = X_sol; % Update guess for continuity along the branch
        
        % Recover the isolated RLE describing function N(jω, A_dc)
        ratio_sol = w * A_delta_c(i) / R_rate;
        N_vals = calc_N(A_delta_c(i), w, R_rate, ratio_sol);
        
        A(i)   = N_vals(1);
        phi(i) = N_vals(2);
    else
        A_delta_c(i) = NaN; phi_sol(i) = NaN; A(i) = NaN; phi(i) = NaN;
    end
end

A_dB    = 20*log10(A);
phi_deg = rad2deg(phi);

% Define logical masks for plotting (using theoretical onset for reference regions)
xi = omega / omega_onset;
m1 = xi < 1;
m2 = (xi >= 1) & (xi < 1.862);
m3 = xi >= 1.862;

% Continuity sanity check (using closest valid solved points)
[~, kI ] = min(abs(xi - 1.000));
[~, kII] = min(abs(xi - 1.862));
fprintf('Boundary check:\n');
fprintf('  ω/ω_onset = 1.000 :  |N| = %+6.3f dB,  ∠N = %+7.3f deg\n', ...
        A_dB(kI),  phi_deg(kI));
fprintf('  ω/ω_onset = 1.862 :  |N| = %+6.3f dB,  ∠N = %+7.3f deg\n\n', ...
        A_dB(kII), phi_deg(kII));

%% =========================================================================
%  FIGURE 1  —  Bode of the DF alone   (≈ Fig. 12)
%% =========================================================================
figure('Color','w','Name','RLE DF — Bode','Position',[80 80 920 620]);

subplot(2,1,1);
semilogx(xi, A_dB, 'b', 'LineWidth', 1.8); grid on; box on;
xline(1.000, '--k', 'LineWidth', 0.9);
xline(1.862, '--k', 'LineWidth', 0.9);
ylabel('|N(j\omega,A)|  [dB]');
title('RLE describing function — harmonic balance solution');

subplot(2,1,2);
semilogx(xi, phi_deg, 'r', 'LineWidth', 1.8); grid on; box on;
xline(1.000, '--k', 'LineWidth', 0.9);
xline(1.862, '--k', 'LineWidth', 0.9);
xlabel('\omega/\omega_{onset}');
ylabel('\angle N(j\omega,A)  [deg]');
ylim([-100 10]);

%% =========================================================================
%  FIGURE 2  —  DF on its own Nichols-style axes (standalone, no matching)
%% =========================================================================
figure('Color','w','Name','RLE DF — Nichols-style','Position',[80 80 760 600]);

hold on;

plot(phi_deg(m1), A_dB(m1), 'b.'); 
plot(phi_deg(m2), A_dB(m2), 'r:');
plot(phi_deg(m3), A_dB(m3), 'go');

plot(0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
text(-3, 1.6, '\omega = \omega_{onset}', 'FontSize', 10);

phi_b = rad2deg(-acos(pi*0.537/2));
A_b   = 20*log10(4*0.537/pi);
plot(phi_b, A_b, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
text(phi_b+1, A_b+1.2, '\omega_{onset}', 'FontSize', 10);
grid on; 
box on;
xlabel('Phase  \angle N(j\omega,A)  [deg]');
ylabel('Gain  |N(j\omega,A)|  [dB]');
title('RLE describing function on Nichols-style axes (standalone)');
legend({'Region I (no sat.)', ...
        'Region II (transition)', ...
        'Region III (full sat.)'}, 'Location','southwest');


%% =========================================================================
%  FIGURE 3  —  Fig. 16 reproduction:
%               Nichols of  G(jω),  N(jω,A),  G·N   with phase matching

% Build all three objects on the SAME omega grid (multiplication well defined)
N_frd  = frd(A .* exp(1j*phi), omega);           % solved describing function alone
G_frd  = frd(Gs_ac_7, omega);                    % linear plant on same grid
LN_frd = G_frd * N_frd;                          % open-loop with RLE DF

%%

figure('Color','w','Name','Fig. 16 — q/q_c with RLE DF','Position',[80 80 900 700]);

hold on;
% 
% nicholsplot(Gs_ac_7, 'b.');

nicholsplot(N_frd, 'r.'); % Piecewise Describing Function's Transfer Function

% nicholsplot(LN_frd, 'k.'); % Aircraft Nonlinear Open-Loop Transfer Function

yline(0,'--'); 

xline(-180,'--'); 

plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r'); % OLOP marker (cursor-picked on the published, matched-phase chart)

text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt), 'Color','r','FontWeight','bold','FontSize',12);

grid on;

legend({'Aircraft''s Linear OLTF','Piecewise Describing Function''s TF','Aircraft''s Nonlinear OLTF','','',''}, 'Location', 'best', 'Orientation', 'horizontal')

xlim auto
ylim auto
% xlim([-300-2 -60-2]);
% ylim([-20-2 15+2]);

hold off;