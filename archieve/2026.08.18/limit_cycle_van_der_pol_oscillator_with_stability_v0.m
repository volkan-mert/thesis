% VAN DER POL OSCILLATOR

% This script:
% 1. solves the Van der Pol oscillator using ode45,
% 2. finds the numerical limit-cycle amplitude and frequency,
% 3. plots the time response and phase portrait,
% 4. plots the describing-function results similar to Slotine Fig. 5.24 and 5.25,
% 5. stores A, w, real(L), and imag(L) in one variable called result.
% 6. checks the stability of the limit cycle by small amplitude disturbances.

clear; clc; close all;

%% Parameters

mu = 1;                 % Van der Pol parameter
tspan = [0 60];         % simulation time
x0 = [0.5 0];           % initial conditions: x(0) = 0.5, xdot(0) = 0

%% Numerical Simulation

% Van der Pol equation:
% xddot + mu*(x^2 - 1)*xdot + x = 0
%
% State variables:
% x(1) = x
% x(2) = xdot

[t,x] = ode45(@(t,x) [x(2); mu*(1-x(1)^2)*x(2)-x(1)], tspan, x0);

% Use the response after 40 seconds as steady state.
i = t >= 40;
t_ss = t(i);
x_ss = x(i,1);
xdot_ss = x(i,2);

%% Numerical Limit-Cycle Amplitude

% For a nearly symmetric oscillation:
% A = (maximum - minimum)/2
A_num = (max(x_ss) - min(x_ss))/2;

%% Numerical Limit-Cycle Frequency

% Find positive-going zero crossings.
z = find(x_ss(1:end-1) < 0 & x_ss(2:end) >= 0);

% Time between consecutive positive-going zero crossings is one period.
T_num = mean(diff(t_ss(z)));

% Convert period to frequency.
f_num = 1/T_num;
w_num = 2*pi*f_num;

% Display numerical results.
fprintf('\nNumerical limit cycle:\n');
fprintf('A = %.4f\n', A_num);
fprintf('w = %.4f rad/s\n', w_num);
fprintf('f = %.4f Hz\n', f_num);
fprintf('T = %.4f s\n\n', T_num);

% Describing-function prediction.
fprintf('Describing function prediction:\n');
fprintf('A = 2\n');
fprintf('w = 1 rad/s\n\n');

%% Time Response

figure(Name='Time Response', NumberTitle='off')
plot(t, x(:,1), 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('x(t)');
title('Time Response', 'Van der Pol Oscillator');

%% Phase Portrait

figure(Name='Phase Portrait', NumberTitle='off')

% Complete trajectory.
plot(x(:,1), x(:,2), 'LineWidth', 1);
hold on;

% Steady-state numerical limit cycle.
plot(x_ss, xdot_ss, 'LineWidth', 2);

% Describing-function prediction: A = 2 and w = 1 rad/s.
numpts = 360;
theta = linspace(0, 2*pi, numpts);

plot(2*sin(theta), 2*cos(theta), 'k--', 'LineWidth', 1.5);

grid on;
axis equal;
xlabel('x');
ylabel('dx/dt');
title('Phase Portrait', 'Van der Pol Oscillator');
legend('Transient trajectory', ...
       'Numerical limit cycle', ...
       'DF prediction: A = 2, \omega = 1', ...
       'Location', 'best');

%% Describing Function

% G(s) = 1/(s^2 - mu*s + 1)
% N(A,w) = j*mu*w*A^2/4
% Limit-cycle condition: G(jw)N(A,w) = -1

% Frequency vector.
w = linspace(0.05, 3, 1000);

% G(jw).
G = 1./((1-w.^2) - 1i*mu*w);

%% Plot Similar to Slotine Fig. 5.24

figure(Name='Slotine''s Book Fig. 5.24', NumberTitle='off')

plot(real(G), imag(G), 'LineWidth', 2);
hold on;

% Amplitude range.
A = linspace(0.5, 4, 400);

% Fixed frequencies.
w0 = [0.5 0.8 1 1.5 2];

for k = 1:length(w0)

    % Calculate N(A,w) for each fixed frequency.
    N = 1i*mu*w0(k)*A.^2/4;

    % Calculate -1/N(A,w).
    M = -1./N;

    % Plot in the complex plane.
    plot(real(M), imag(M), '--', 'LineWidth', 1);

end

% Mark the theoretical solution A = 2, w = 1 rad/s.
plot(0, 1/mu, 'ko', 'MarkerFaceColor', 'k');
text(0.05, 1/mu + 0.5, 'A = 2, \omega = 1');

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Nyquist Plot');

txt1 = {'Van der Pol Oscillator Limit Cycle Detection for different frequencies', ...
        '(Slotine''s Book Fig. 5.24)'};
subtitle(txt1);

legend('G(j\omega), Linear TF', ...
       '-1/N(A,\omega), \omega = 0.5', ...
       '-1/N(A,\omega), \omega = 0.8', ...
       '-1/N(A,\omega), \omega = 1', ...
       '-1/N(A,\omega), \omega = 1.5', ...
       '-1/N(A,\omega), \omega = 2', ...
       'Limit Cycle', ...
       'Location', 'best');

%% Plot Similar to Slotine Fig. 5.25

% This plot shows G(jw)N(A,w) for different fixed amplitudes.
% A limit cycle occurs when the curve reaches (-1,0).

figure(Name='Slotine''s Book Fig. 5.25', NumberTitle='off')
hold on;

% Fixed amplitude values.
A_values = [1 1.5 2 2.5 3];

% result columns:
% Column 1 = A
% Column 2 = w
% Column 3 = real(L)
% Column 4 = imag(L)
result = [];

for k = 1:length(A_values)

    % Current amplitude.
    A = A_values(k);

    % Calculate N(A,w).
    N = 1i*mu*w*A^2/4;

    % Calculate L = G(jw)N(A,w).
    L = G.*N;

    % Store A, w, real(L), and imag(L) in one variable.
    result = [result;
              A*ones(length(w),1), w(:), real(L(:)), imag(L(:))];

    % Plot the result.
    plot(real(L), imag(L), 'LineWidth', 1.5);

end

%% Find the Point Closest to (-1,0)

% Calculate the distance of every point from (-1,0).
distance = sqrt((result(:,3) + 1).^2 + result(:,4).^2);

% Find the smallest distance.
[minimum_distance, index] = min(distance);

% Get A, w, real(L), and imag(L) at this point.
A_limit = result(index,1);
w_limit = result(index,2);
real_limit = result(index,3);
imag_limit = result(index,4);

% Display the result.
fprintf('Point closest to (-1,0):\n');
fprintf('A = %.4f\n', A_limit);
fprintf('w = %.4f rad/s\n', w_limit);
fprintf('real(L) = %.6f\n', real_limit);
fprintf('imag(L) = %.6f\n', imag_limit);
fprintf('distance from (-1,0) = %.6f\n\n', minimum_distance);

%% Complete Fig. 5.25 Plot

% Mark the critical point (-1,0).
plot(-1, 0, 'ko', 'MarkerFaceColor', 'k');
text(-0.95, 0.05, '(-1,0): Limit Cycle');

% Mark the closest calculated point.
plot(real_limit, imag_limit, 'ro', 'MarkerFaceColor', 'r');
text(real_limit + 0.05, imag_limit - 0.05, ...
    sprintf('A = %.2f, \omega = %.4f', A_limit, w_limit));

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Nyquist Plot');

txt2 = {'Van der Pol Oscillator Limit Cycle Detection for different amplitudes', ...
        '(Slotine''s Book Fig. 5.25)'};
subtitle(txt2);

legend('A = 1', ...
       'A = 1.5', ...
       'A = 2', ...
       'A = 2.5', ...
       'A = 3', ...
       '(-1,0): Limit Cycle', ...
       'Closest calculated point', ...
       'Location', 'best');

%% Stability of the Limit Cycle

% Slotine Fig. 5.26 explains stability by slightly increasing or decreasing
% the limit-cycle amplitude. A stable limit cycle should return to the
% original oscillation after a small disturbance.
%
% The simple "encircled / not encircled" explanation in Fig. 5.26 assumes
% that G(s) has no unstable poles. In this Van der Pol example,
% G(s) = 1/(s^2 - mu*s + 1) has unstable poles. Therefore, the disturbance
% behavior is checked directly here by numerical simulation.

% Small change in amplitude around the predicted limit cycle.
dA = 0.25;

A_low = A_limit - dA;
A_high = A_limit + dA;

% Initial conditions below and above the limit-cycle amplitude.
x0_low = [A_low 0];
x0_high = [A_high 0];

% Run two additional simulations.
[t_low,x_low] = ode45(@(t,x) [x(2); mu*(1-x(1)^2)*x(2)-x(1)], tspan, x0_low);

[t_high,x_high] = ode45(@(t,x) [x(2); mu*(1-x(1)^2)*x(2)-x(1)], tspan, x0_high);

% Use the last 20 seconds to calculate the final amplitudes.
i_low = t_low >= 40;
i_high = t_high >= 40;

A_end_low = (max(x_low(i_low,1)) - min(x_low(i_low,1)))/2;
A_end_high = (max(x_high(i_high,1)) - min(x_high(i_high,1)))/2;

% Compare the final amplitudes with the numerical limit-cycle amplitude.
tolerance = 0.05*A_num;

if abs(A_end_low - A_num) < tolerance && abs(A_end_high - A_num) < tolerance

    stability = 'STABLE';

else

    stability = 'UNSTABLE';

end

% Display the stability result.
fprintf('Limit-cycle stability check:\n');
fprintf('Initial amplitude below limit cycle = %.4f\n', A_low);
fprintf('Final amplitude from lower disturbance = %.4f\n', A_end_low);
fprintf('Initial amplitude above limit cycle = %.4f\n', A_high);
fprintf('Final amplitude from upper disturbance = %.4f\n', A_end_high);
fprintf('The limit cycle is %s.\n\n', stability);


%% Stability Plot Similar to Slotine Fig. 5.26

figure(Name='Slotine''s Book Fig. 5.26', NumberTitle='off')

% Plot G(jw).
plot(real(G), imag(G), 'LineWidth', 2);
hold on;

% Plot -1/N(A,w) near the limit-cycle frequency.
A_stability = linspace(0.8*A_limit, 1.2*A_limit, 300);
N_stability = 1i*mu*w_limit*A_stability.^2/4;
M_stability = -1./N_stability;

plot(real(M_stability), imag(M_stability), '--', 'LineWidth', 1.5);

% Calculate the three points: smaller A, limit-cycle A, and larger A.
M_low = -1/(1i*mu*w_limit*A_low^2/4);
M_limit = -1/(1i*mu*w_limit*A_limit^2/4);
M_high = -1/(1i*mu*w_limit*A_high^2/4);

plot(real(M_low), imag(M_low), 'bo', 'MarkerFaceColor', 'b');
plot(real(M_limit), imag(M_limit), 'ko', 'MarkerFaceColor', 'k');
plot(real(M_high), imag(M_high), 'ro', 'MarkerFaceColor', 'r');

text(real(M_low)+0.03, imag(M_low), 'A decreased');
text(real(M_limit)+0.03, imag(M_limit), 'Limit cycle');
text(real(M_high)+0.03, imag(M_high), 'A increased');

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Stability of Limit Cycle');
subtitle('Small amplitude disturbances around the limit-cycle point');

legend('G(j\omega)', ...
       '-1/N(A,\omega)', ...
       'A < A_{limit}', ...
       'A = A_{limit}', ...
       'A > A_{limit}', ...
       'Location', 'best');


%% Phase-Plane Verification of Stability

figure(Name='Limit Cycle Stability Verification', NumberTitle='off')

% Trajectory starting below the limit-cycle amplitude.
plot(x_low(:,1), x_low(:,2), 'LineWidth', 1);
hold on;

% Trajectory starting above the limit-cycle amplitude.
plot(x_high(:,1), x_high(:,2), 'LineWidth', 1);

% Numerical steady-state limit cycle.
plot(x_ss, xdot_ss, 'k', 'LineWidth', 2);

grid on;
axis equal;
xlabel('x');
ylabel('dx/dt');
title('Limit Cycle Stability Verification');
subtitle(['The limit cycle is ', stability]);

legend('Initial amplitude below limit cycle', ...
       'Initial amplitude above limit cycle', ...
       'Steady-state limit cycle', ...
       'Location', 'best');

