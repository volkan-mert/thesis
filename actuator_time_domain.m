% Rate Limiting Element (Hysteresis and Saturation)

clear; clc; close all;

% --- 1. Actuator Parameters ---
K = 20;   % Forward gain
S = 15;   % Upper saturation limit (max rate)
R = -S;  % Lower saturation limit (min rate)

% --- 2. Simulation Setup ---
tspan = [0 5]; % Simulate for 5 seconds
x0 = 0;        % Initial actuator position

% Pilot input command (Sine wave with amplitude 2 to force saturation)
pilot_cmd = @(t) 2 * sin(pi * t); 

% --- 3. Equations of Motion ---

actuator_dyn = @(t, x) max(R, min(S, K * (pilot_cmd(t) - x))); % The min/max logic naturally clips the rate to S and R.

% --- 4. Numerical Integration ---
[t, x] = ode45(actuator_dyn, tspan, x0);

% --- 5. Post-Processing ---
% Recalculate input and rate over the solved time steps for plotting
u = pilot_cmd(t);
x_dot = max(R, min(S, K .* (u - x))); 

% --- 6. Plotting ---
figure('Name', 'Actuator Response');

% Top Plot: Position Tracking
subplot(2,1,1);
plot(t, u, '--k', 'LineWidth', 1.5); hold on;
plot(t, x, 'b', 'LineWidth', 1.5);
grid on;
title('Actuator Position Tracking');
xlabel('Time (sec)');
ylabel('Position');
legend('Pilot Command', 'Actuator State', 'Location', 'best');

% Bottom Plot: Rate Saturation
subplot(2,1,2);
plot(t, x_dot, 'r', 'LineWidth', 1.5); hold on;
yline(S, '--k', 'Upper Limit');
yline(R, '--k', 'Lower Limit');
grid on;
title('Actuator Rate (Showing Saturation Limits)');
xlabel('Time (sec)');
ylabel('Rate');
ylim([-20 20]);