% VAN DER POL OSCILLATOR
%
% This script:
% 1. solves the Van der Pol oscillator with ode45,
% 2. estimates the numerical limit-cycle amplitude and frequency,
% 3. plots the time response and phase portrait,
% 4. draws describing-function plots similar to Slotine Fig. 5.24 and 5.25,
% 5. stores result = [A, w, real(L), imag(L)] for Fig. 5.25,
% 6. finds the point closest to (-1,0).

clear; clc; close all;

%% Parameters for Numerical Simulation

% Van der Pol parameter
mu = 1;

% Simulation time
tspan = [0 60];

% Initial conditions
x0 = [0.5 0];   % x(0)=0.5, xdot(0)=0

%% Numerical simulation of the Van der Pol Oscillator

% Van der Pol equation:
% xddot + mu*(x^2 - 1)*xdot + x = 0
%
% State variables:
% x(1) = x
% x(2) = dx/dt

[t,x] = ode45(@(t,x) [x(2); mu*(1-x(1)^2)*x(2)-x(1)], tspan, x0);

% Use steady-state part after 40 seconds
i = t >= 40;
t_ss    = t(i);
x_ss    = x(i,1);
xdot_ss = x(i,2);

% Numerical amplitude
A_num = (max(x_ss)-min(x_ss))/2;

% Numerical frequency: positive-going zero crossings
z = find(x_ss(1:end-1) < 0 & x_ss(2:end) >= 0);

% Numerical period and frequency
T_num = mean(diff(t_ss(z)));
f_num = 1/T_num;
w_num = 2*pi*f_num;

% Display numerical limit-cycle results
fprintf('\nNumerical limit cycle:\n');
fprintf('A = %.4f\n', A_num);
fprintf('w = %.4f rad/s\n', w_num);
fprintf('f = %.4f Hz\n', f_num);
fprintf('T = %.4f s\n\n', T_num);

% Describing-function prediction
fprintf('Describing function prediction:\n');
fprintf('A = 2\n');
fprintf('w = 1 rad/s\n\n');

%% Time response

figure('Name','Time Response','NumberTitle','off');
plot(t, x(:,1), 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('x(t)');
title('Time Response');
subtitle('Van der Pol Oscillator');

%% Phase portrait

figure('Name','Phase Portrait','NumberTitle','off');
plot(x(:,1), x(:,2), 'LineWidth', 1);
hold on;
plot(x_ss, xdot_ss, 'LineWidth', 2);

% Describing-function prediction: A = 2, w = 1
numpts = 360;
theta = linspace(0, 2*pi, numpts);
plot(2*sin(theta), 2*cos(theta), 'k--', 'LineWidth', 1.5);

grid on;
axis equal;
xlabel('x');
ylabel('dx/dt');
title('Phase Portrait');
subtitle('Van der Pol Oscillator');
legend('Transient Trajectory', ...
       'Steady-State Limit Cycle', ...
       'Theoretical Limit Cycle (DF Prediction)', ...
       'Location','best');

%% Describing function generation
% G(s) = 1 / (s^2 - mu*s + 1)
% N(A,w) = j*mu*w*A^2/4
% Limit cycle when G(jw)N(A,w) = -1

% Frequency vector
w = linspace(0.05, 3, 1000);

% Linear frequency response G(jw)
G = 1 ./ ((1 - w.^2) - 1i*mu*w);

%% Plot similar to Slotine Fig. 5.24
% Plot G(jw) and a family of -1/N(A,w) curves for different fixed w

figure('Name','Slotine Fig. 5.24','NumberTitle','off');
plot(real(G), imag(G), 'LineWidth', 2);
hold on;

% Amplitude range
A = linspace(0.5, 4, 400);

% Fixed frequency values
w0 = [0.5 0.8 1 1.5 2];

for k = 1:length(w0)
    % Describing function for current fixed frequency
    N = 1i*mu*w0(k)*A.^2/4;

    % Negative inverse describing function
    M = -1 ./ N;

    % Plot in complex plane
    plot(real(M), imag(M), '--', 'LineWidth', 1);
end

% Mark the DF solution A = 2, w = 1
plot(0, 1/mu, 'ko', 'MarkerFaceColor', 'k');
text(0.05, 1/mu + 0.1, 'A = 2, \omega = 1');

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Nyquist Plot');
subtitle({'Van der Pol Oscillator Limit Cycle Detection for Different Frequencies', ...
          '(Similar to Slotine Fig. 5.24)'});
legend('G(j\omega), Linear TF', ...
       '-1/N(A,\omega) at \omega = 0.5', ...
       '-1/N(A,\omega) at \omega = 0.8', ...
       '-1/N(A,\omega) at \omega = 1', ...
       '-1/N(A,\omega) at \omega = 1.5', ...
       '-1/N(A,\omega) at \omega = 2', ...
       'Limit Cycle Point', ...
       'Location','best');

%% Plot similar to Slotine Fig. 5.25
% Plot G(jw)N(A,w) for several fixed amplitude values
% result = [A, w, real(L), imag(L)]

figure('Name','Slotine Fig. 5.25','NumberTitle','off');
hold on;

% Fixed amplitude values
A_values = [1 1.5 2 2.5 3];

% Single storage variable:
% Columns = [A, w, real(L), imag(L)]
result = [];

for k = 1:length(A_values)

    % Current amplitude
    A = A_values(k);

    % Describing function
    N = 1i*mu*w*A^2/4;

    % Open-loop complex value
    L = G .* N;

    % Store [A, w, real(L), imag(L)]
    result = [result;
              A*ones(length(w),1), w(:), real(L(:)), imag(L(:))];

    % Plot the curve
    plot(real(L), imag(L), 'LineWidth', 1.5);
end

% Find the point closest to (-1,0)
% For finding the point closest to (a,b), distance = sqrt((Real{L} - a).^2 + (Imaginary{L} - b).^2)
distance = sqrt((result(:,3) + 1).^2 + result(:,4).^2);
[minimum_distance, index] = min(distance);

A_limit    = result(index,1);
w_limit    = result(index,2);
real_limit = result(index,3);
imag_limit = result(index,4);

% Display result in command window
fprintf('Point closest to (-1,0):\n');
fprintf('A = %.4f\n', A_limit);
fprintf('w = %.4f rad/s\n', w_limit);
fprintf('real(L) = %.6f\n', real_limit);
fprintf('imag(L) = %.6f\n', imag_limit);
fprintf('distance from (-1,0) = %.6e\n\n', minimum_distance);

% Mark the critical point (-1,0)
plot(-1, 0, 'ko', 'MarkerFaceColor', 'k');
text(-0.95, 0.05, '(-1,0): Limit Cycle');

% Mark the closest calculated point
plot(real_limit, imag_limit, 'ro', 'MarkerFaceColor', 'r');
text(real_limit + 0.05, imag_limit - 0.05, ...
    sprintf('A = %.2f, \\omega = %.4f', A_limit, w_limit));

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Nyquist Plot');
subtitle({'Van der Pol Oscillator Limit Cycle Detection for Different Amplitudes', ...
          '(Similar to Slotine Fig. 5.25)'});
legend('A = 1', ...
       'A = 1.5', ...
       'A = 2', ...
       'A = 2.5', ...
       'A = 3', ...
       '(-1,0): Limit Cycle Crossing', ...
       'Closest Calculated Point', ...
       'Location','best');

%% Notes
% result(:,1) = A
% result(:,2) = w
% result(:,3) = real(L)
% result(:,4) = imag(L)
%
% Example:
% result(1:10,:)
% will show the first 10 rows of [A, w, real(L), imag(L)].