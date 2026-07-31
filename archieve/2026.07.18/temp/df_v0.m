%% Negative Inverse Describing Function Analysis & Simulation
% Predicts and verifies limit cycles for a linear plant with relay feedback.
clear; clc; close all;

%% 1. System & Nonlinearity Definition
K = 10;                     % Linear plant gain
M = 1;                      % Relay output amplitude (Ideal Relay)

% Linear Plant: G(s) = K / (s * (s + 1) * (s + 2))
% Characteristic denominator polynomial: s^3 + 3*s^2 + 2*s
num = K;
den = [1, 3, 2, 0];

%% 2. Frequency-Domain Analysis (Describing Function Technique)
% Generate frequency vector
w = logspace(-1, 1.5, 1000);
s = 1i * w;

% Evaluate plant frequency response G(jw)
G_jw = K ./ (s .* (s + 1) .* (s + 2));

% Describing Function for Ideal Relay: N(A) = 4*M / (pi*A)
% Negative Inverse: -1/N(A) = -pi*A / (4*M)
A = linspace(0.1, 5, 500);          % Sinusoidal input amplitude vector
neg_inv_N = -(pi * A) / (4 * M);    % Trajectory lies along the negative real axis

% Find Analytical Intersection (Phase Crossover)
% Im{G(jw)} = 0 when 2 - w^2 = 0 -> w = sqrt(2)
w_pc = sqrt(2);
G_pc = K / (-3 * w_pc^2);           % Value of G(jw) at phase crossover (-K/6)

% Calculate predicted limit cycle amplitude from intersection: G(jw) = -1/N(A)
A_pred = abs(G_pc) * (4 * M) / pi;

fprintf('--- Describing Function Prediction ---\n');
fprintf('Predicted Limit Cycle Frequency: %.4f rad/s (%.4f Hz)\n', w_pc, w_pc/(2*pi));
fprintf('Predicted Limit Cycle Amplitude: %.4f\n\n', A_pred);

%% 3. Time-Domain Verification (ode45)
% State-space realization: x1 = y, x2 = y_dot, x3 = y_ddot
A_sys = [0,  1,  0;
    0,  0,  1;
    0, -2, -3];
B_sys = [0; 0; K];
C_sys = [1, 0, 0];

% Closed-loop dynamics with feedback relay: u(t) = -M * sign(y(t))
ode_func = @(t, x) A_sys * x + B_sys * (-M * sign(C_sys * x));

% Simulate time response
t_span = [0, 25];
x0 = [0.5; 0; 0];                   % Small initial disturbance to trigger oscillation
[t, X] = ode45(ode_func, t_span, x0);
y = X(:, 1);

% Extract steady-state limit cycle metrics (without Signal Processing Toolbox)
idx_ss = find(t > 15);
y_ss = y(idx_ss);
t_ss = t(idx_ss);

% Locate local peaks in steady state
dy = diff(y_ss);
peak_idx = find(dy(1:end-1) > 0 & dy(2:end) <= 0) + 1;
A_sim = mean(y_ss(peak_idx));
w_sim = 2 * pi / mean(diff(t_ss(peak_idx)));

fprintf('--- Time-Domain Simulation Results ---\n');
fprintf('Simulated Limit Cycle Frequency: %.4f rad/s\n', w_sim);
fprintf('Simulated Limit Cycle Amplitude: %.4f\n', A_sim);

%% 4. Graphical Visualization
figure('Name', 'Describing Function Limit Cycle Analysis', 'Color', 'w', 'Position', [100, 100, 900, 400]);

% Nyquist Plot & -1/N(A) Locus
subplot(1, 2, 1);
plot(real(G_jw), imag(G_jw), 'b-', 'LineWidth', 1.5); hold on;
plot(real(G_jw), -imag(G_jw), 'b--', 'LineWidth', 1.0); % Mirror for negative freq
plot(neg_inv_N, zeros(size(neg_inv_N)), 'r-', 'LineWidth', 2);
plot(G_pc, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Intersection point
grid on;
axis([-2.5, 0.5, -2, 2]);
xlabel('Real Axis'); ylabel('Imaginary Axis');
title('Nyquist Plot vs. -1/N(A) Locus');
legend('G(j\omega)', 'G(-j\omega)', '-1/N(A) Locus', 'Limit Cycle Intersection', 'Location', 'northwest');

% Time-Domain Response
subplot(1, 2, 2);
plot(t, y, 'k-', 'LineWidth', 1.2); hold on;
yline(A_pred, 'r--', 'Predicted Bound', 'LineWidth', 1.2);
yline(-A_pred, 'r--', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('Output y(t)');
title('Time-Domain Verification (Limit Cycle)');
legend('Simulated y(t)', 'DF Predicted Bound', 'Location', 'southeast');