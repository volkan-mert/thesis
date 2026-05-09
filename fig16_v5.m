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
[w_opt, ~] = fmincon(f, 5, [],[],[],[], 0, 100, [], fminopt);

fprintf('Closed-loop onset frequency  ω_onset = %.4f rad/s\n\n', w_opt);

%% =========================================================================
%  Piecewise RLE describing function   (Duda / Hanke / Gilbreath)
%% =========================================================================
omega_onset = w_opt;                             % from Eq. 18  — NOT 1
% omega       = logspace(-1, 2, 1000);              % must cover ω = 1 rad/s
omega = 0.1:0.1:10;
xi          = omega / omega_onset;               % ω / ω_onset

A   = zeros(size(omega));                        % |N(jω,A)|
phi = zeros(size(omega));                        % ∠N(jω,A)  [rad]

% Define logical masks for plotting later
m1 = xi < 1;
m2 = (xi >= 1) & (xi < 1.862);
m3 = xi >= 1.862;

% Calculate describing function element by element
for i = 1:size(omega, 2)
    if xi(i) < 1
        % Region I  —  no saturation
        A(i)   = 1;
        phi(i) = 0;

        A1(i) = A(i);
        phi1(i) = i;

    elseif xi(i) >= 1 && xi(i) < 1.862
        % Region II —  transition (cubic in α = ω/ω_onset)
        a = xi(i);
        A(i)   = 0.2908*a^3 - 1.4396*a^2 + 1.9232*a + 0.2230;
        phi(i) = 0.5280*a^3 - 2.6213*a^2 + 3.5056*a - 1.4171;

        A2(i) = A(i);
        phi2(i) = i;

    elseif xi(i) >= 1.862
        % Region III —  fully saturated  (ϖ = ω_onset/ω)
        v = omega_onset / omega(i);
        A(i)   = 4*v/pi;
        phi(i) = -acos(pi*v/2);

        A3(i) = A(i);
        phi3(i) = i;

    end
end



% plot(xi); hold; plot(omega)
hold on
% plot(x1*0.1, r1)
% hold on
% plot(x2*0.1, r2)
% hold on
% plot(x3*0.1, r3)
% hold on
% xline((w_opt-0.1)/(10-0.1).*10,'Color','r','LineStyle','-.')
legend show
% 
%         mask2 = (x2 ~= 0);
%         mask3 = (x3 ~= 0);
%         mask1 = (x1 ~= 0);
% 
% plot(x1(mask1)*0.1, r1(mask1), 'r--o', x2(mask2)*0.1, r2(mask2), 'k-->', x3(mask3)*0.1, r3(mask3), 'c--')
% 
% max(r1)
% max(r2)
% max(r3)

% 
% figure;
% hold on
% plot(A1);
% plot(A2);
% plot(A3);
% hold off
% legend show
% 
% figure;
% hold on
% plot(phi1);
% plot(phi2);
% plot(phi3);
% hold off
% legend show

% sys_frd1 = frd(A1, phi1);
% sys_frd2 = frd(A2, phi2);
% sys_frd3 = frd(A3, phi3);

opts = nicholsoptions;
opts.PhaseMatching       = 'on';
opts.PhaseMatchingFreq   = 1;            % rad/s
opts.PhaseMatchingValue  = -180;         % anchor phase at -180° at ω = 1 rad/s
opts.PhaseWrapping       = 'on';
opts.PhaseWrappingBranch = -360;

sys_frd1 = frd(A1, phi1);
sys_frd2 = frd(A2, phi2);
sys_frd3 = frd(A3, phi3);
% 
% % Build all three objects on the SAME omega grid (multiplication well defined)
% N_frd1  = frd(A1 .* exp(1j*phi1), phi1);           % describing function alone
% N_frd2  = frd(A2 .* exp(1j*phi2), phi2);           % describing function alone
% N_frd3  = frd(A3 .* exp(1j*phi3), phi3);           % describing function alone
% 
% G_frd  = frd(Gs_ac_7, omega);                    % linear plant on same grid
% 
% LN_frd1 = G_frd * N_frd1;                          % open-loop with RLE DF
% LN_frd2 = G_frd * N_frd2;                          % open-loop with RLE DF
% LN_frd3 = G_frd * N_frd3;                          % open-loop with RLE DF
% 
% figure('Color','w','Name','Fig. 16 — q/q_c with RLE DF','Position',[80 80 900 700]);
% hold on;
% nicholsplot(Gs_ac_7, opts);   
% nicholsplot(N_frd1,   opts);
% nicholsplot(N_frd2,   opts);
% nicholsplot(N_frd3,   opts);
% nicholsplot(LN_frd1,  opts);
% nicholsplot(LN_frd2,  opts);
% nicholsplot(LN_frd3,  opts);
% %
% yline(0,'--');
% xline(-180,'--');
% 
% % OLOP marker (cursor-picked on the published, matched-phase chart)
% plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop, 'rp', ...
%      'MarkerSize', 15, 'MarkerFaceColor', 'r');
% text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, ...
%      sprintf('\\omega_{OLOP} = %.3f rad/s', w_opt), ...
%      'Color','r','FontWeight','bold','FontSize',12);
% 
% grid on;
% xlim([-300-2 -60-2]);
% ylim([-20-2 15+2]);
% legend({'G_{ac\_7}', ...
%         'N(j\omega,A)', ...
%         'G_{ac\_7}\cdot N(j\omega,A)', ...
%         'OLOP'}, 'Location','southwest');
% hold off;
