%% ================================================================
%%  Example 4.6-1  —  Simple Pitch-Attitude Hold Autopilot
%%  Transport Aircraft: 25 000 ft, 500 ft/s TAS, xcg = 0.25c
%%  Reference: Stevens & Lewis, "Aircraft Control and Simulation"
%%
%%  LOOP STRUCTURE (Figure 4.6-1):
%%
%%   theta_c -->[ ]--> Gc(s) -->[ ]--> Actuator --> A/C --> theta
%%              ^ e      (=1)    ^ u1   -10/(s+10)         |
%%              |                |                           |
%%              |<----  k_th <---|<---  kq <-- Rate gyro <--|--- q
%%              |                                           |
%%              |<------- Attitude gyro <------------------|
%%
%%  Open-Loop  : L(s)    = k_th * [theta/u1]_inner_CL
%%  Closed-Loop: theta/theta_c = L(s) / (1 + L(s))
%% ================================================================
clear; clc; close all;

%% ---------------------------------------------------------------
%%  1.  Aircraft Plant  (5 states: vT, alpha, theta, q, h)
%% ---------------------------------------------------------------
A = [-0.0082354   18.938   -32.170   0.0        5.9022e-5 ;
     -0.00025617  -0.56761   0.0      1.0        2.2633e-6 ;
      0.0          0.0       0.0      1.0        0.0       ;
      1.3114e-5   -1.4847    0.0     -0.47599   -1.4947e-7 ;
      0.0        -500.0    500.0      0.0        0.0       ];

B = [0; 0; 0; -0.019781; 0];    % single input: delta_e  [rad]

%  Outputs (rad->deg conversion via 57.296)
%  Row 1: theta [deg]
%  Row 2: q     [deg/s]
C = [0  0  57.296   0       0 ;   % theta
     0  0   0      57.296   0 ];  % q
D = zeros(2,1);

plant = ss(A, B, C, D);

%% ---------------------------------------------------------------
%%  2.  Elevator Actuator  (simple lag, tau = 0.1 s)
%%      TF = -10/(s+10)   [sign change: +u -> +delta_e -> +pitchrate]
%% ---------------------------------------------------------------
actua = ss(-10, 10, -1, 0);   % ca = -1 handles sign convention

%% ---------------------------------------------------------------
%%  3.  Combined system: Actuator + Aircraft
%%      Input  : u_1  (control signal before actuator)
%%      Outputs: [theta ; q]   (both in degrees)
%% ---------------------------------------------------------------
sys1 = series(actua, plant);
[ap, bp, cp, dp] = ssdata(sys1);
%  cp(1,:) → theta output   (output index 1)
%  cp(2,:) → q     output   (output index 2)

%% ---------------------------------------------------------------
%%  4.  Design Gains (from Example 4.6-1)
%% ---------------------------------------------------------------
kq   = 2.5;    % inner rate  feedback  [elevator-deg / (deg/s pitch rate)]
k_th = 4.0;    % outer attitude feedback [elevator-deg / deg pitch]

%% ---------------------------------------------------------------
%%  5.  Close Inner Rate Loop  (q feedback with gain kq)
%%      u_1 = u_outer - kq * q
%%      => A_innerCL = ap - bp * kq * cp(2,:)
%% ---------------------------------------------------------------
A_icl = ap - bp * kq * cp(2,:);      % inner closed-loop A matrix

%  Inner CL TF:  theta / u_1   (transfer function seen by outer loop)
G_in_CL = ss(A_icl, bp, cp(1,:), dp(1));
G_in_CL = minreal(G_in_CL, 1e-3);   % cancel near pole-zero pairs

%% ---------------------------------------------------------------
%%  6.  Open-Loop Transfer Function  L(s)
%%      Breaking the outer attitude loop at theta_c:
%%        L(s) = k_th * G_in_CL(s)
%% ---------------------------------------------------------------
L = k_th * G_in_CL;
L = minreal(L, 1e-3);

%% ---------------------------------------------------------------
%%  7.  Closed-Loop Transfer Function  theta / theta_c
%%      Unity negative feedback on the outer attitude loop
%% ---------------------------------------------------------------
G_CL = feedback(L, 1);

%% ---------------------------------------------------------------
%%  8.  Book Closed-Loop TF — Equation (3)
%%      Poles:  -1.999 ± j2.389  (short period, damped)
%%              -6.646           (actuator-related)
%%              -0.3815
%%              -0.02522         (phugoid, damped)
%%              -1.718e-4        (altitude, near-cancel)
%%      Zeros:  -0.5567,  -0.01897,  -1.666e-4
%% ---------------------------------------------------------------
sp_quad = [1,  2*1.999,  1.999^2 + 2.389^2];   % (s+1.999)^2 + 2.389^2
num3 = 45.33 * conv([1, 0.5567], conv([1, 0.01897], [1, 1.666e-4]));
den3 = conv( conv(sp_quad, [1, 6.646]), ...
             conv([1, 0.3815], conv([1, 0.02522], [1, 1.718e-4])) );
G_CL_book = tf(num3, den3);

%% ---------------------------------------------------------------
%%  9.  Display TFs and Stability Margins
%% ---------------------------------------------------------------
fprintf('\n===========================================================\n');
fprintf('  TRANSFER FUNCTIONS — Example 4.6-1 Pitch-Attitude Hold\n');
fprintf('===========================================================\n');

fprintf('\n--- Combined Plant+Actuator TF :  theta / delta_e ---\n');
disp( zpk(minreal(ss(ap, bp, cp(1,:), dp(1)), 1e-3)) );

fprintf('\n--- Inner CL TF :  theta / u_1  (kq = %.1f) ---\n', kq);
disp( zpk(G_in_CL) );

fprintf('\n--- Open-Loop TF :  L(s) = k_th * (theta/u1)_innerCL ---\n');
fprintf('    k_th = %.1f,  kq = %.1f\n', k_th, kq);
disp( zpk(L) );

fprintf('\n--- Closed-Loop TF :  theta / theta_c  [Computed] ---\n');
disp( zpk(G_CL) );

fprintf('\n--- Closed-Loop TF :  theta / theta_c  [Book Eq.(3)] ---\n');
disp( zpk(G_CL_book) );

[Gm, Pm, Wgm, Wpm] = margin(L);
fprintf('\n--- Stability Margins (Open Loop) ---\n');
fprintf('  Gain  Margin = %7.2f dB   at  omega = %.4f rad/s\n', 20*log10(Gm), Wgm);
fprintf('  Phase Margin = %7.2f deg  at  omega = %.4f rad/s\n', Pm, Wpm);
fprintf('  DC Gain of CL = %.4f  (Book value ≈ 0.77)\n', dcgain(G_CL));

%% ---------------------------------------------------------------
%%  PLOT 1:  Root Locus — Rate Gain kq  (Figure 4.6-2)
%%           Outer loop closed at k_th = 4.0, sweeping kq
%% ---------------------------------------------------------------
figure('Name','Root Locus — kq (inner rate feedback)','Color','w','Position',[50 550 600 500]);

A_outer = ap - bp * k_th * cp(1,:);     % outer loop closed at k_th=4, inner open
sys_rl  = ss(A_outer, bp, cp(2,:), dp(2));   % open inner: q output

rlocus(sys_rl);
hold on;
% Mark chosen gain kq = 2.5
cl_poles_kq = eig(A_outer - bp * kq * cp(2,:));
plot(real(cl_poles_kq), imag(cl_poles_kq), 'gs', 'MarkerSize', 10, ...
     'MarkerFaceColor','g', 'DisplayName', sprintf('k_q = %.1f', kq));
hold off;
xlim([-12 2]); ylim([-6 6]);
title(sprintf('Root Locus — Rate Feedback k_q  (k_\\theta = %.1f fixed)', k_th));
xlabel('Real Axis'); ylabel('Imaginary Axis');
legend('show','Location','northeast');
grid on;

%% ---------------------------------------------------------------
%%  PLOT 2:  Bode Diagram
%% ---------------------------------------------------------------
figure('Name','Bode Diagram — Open Loop L(s)','Color','w','Position',[700 550 700 550]);

bodeOpts = bodeoptions;
bodeOpts.FreqUnits = 'rad/s';
bodeOpts.Grid      = 'on';
bode(L, {1e-3, 100}, bodeOpts);

% Overlay stability margins
margin(L);    % adds Gm/Pm annotations automatically
title(sprintf('Bode Diagram — L(s),   k_q = %.1f,  k_\\theta = %.1f', kq, k_th));

%% ---------------------------------------------------------------
%%  PLOT 3:  Nyquist Plot
%% ---------------------------------------------------------------
figure('Name','Nyquist Plot — Open Loop L(s)','Color','w','Position',[50 50 600 500]);

nyquist(L);
grid on;
title('Nyquist Plot — Open Loop L(s)');
% Mark critical point
hold on;
plot(-1, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Critical Point (-1,0)');
hold off;
legend('show','Location','best');

%% ---------------------------------------------------------------
%%  PLOT 4:  Nichols Chart
%% ---------------------------------------------------------------
figure('Name','Nichols Chart — Open Loop L(s)','Color','w','Position',[700 50 700 550]);

nichols(L, {1e-3, 100});
ngrid;
title(sprintf('Nichols Chart — L(s),   k_q = %.1f,  k_\\theta = %.1f', kq, k_th));
grid on;

%% ---------------------------------------------------------------
%%  PLOT 5:  Step Response — Closed Loop  theta/theta_c
%% ---------------------------------------------------------------
figure('Name','Step Response — Closed Loop','Color','w','Position',[400 300 650 450]);

t = 0:0.1:50;
step(G_CL,      t); hold on;
step(G_CL_book, t, 'r--');
yline(dcgain(G_CL), 'k:', 'LineWidth',1);
hold off;
legend('Computed (state-space)', 'Book Eq.(3)', ...
       sprintf('SS = %.3f', dcgain(G_CL)), 'Location','best');
xlabel('Time (s)');
ylabel('\theta / \theta_c  (amplitude)');
title('Step Response — Closed Loop   \theta / \theta_c');
grid on;

%% ---------------------------------------------------------------
%%  PLOT 6:  All three frequency plots on one figure (summary)
%% ---------------------------------------------------------------
figure('Name','Frequency Response Summary','Color','w','Position',[200 100 1100 400]);

subplot(1,3,1);
bode(L, {1e-3,100}); grid on;
title('Bode');

subplot(1,3,2);
nyquist(L); grid on;
hold on; plot(-1,0,'rx','MarkerSize',10,'LineWidth',2); hold off;
title('Nyquist');

subplot(1,3,3);
nichols(L, {1e-3,100}); ngrid; grid on;
title('Nichols');

sgtitle(sprintf('Open-Loop L(s)  |  k_q = %.1f,  k_\\theta = %.1f', kq, k_th), ...
        'FontSize', 13, 'FontWeight','bold');
