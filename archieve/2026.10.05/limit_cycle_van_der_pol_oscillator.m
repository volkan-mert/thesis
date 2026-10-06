% VAN DER POL OSCILLATOR

% This script:
% 1. solves the Van der Pol oscillator with ode45,
% 2. estimates the numerical limit-cycle amplitude and frequency,
% 3. plots the time response and phase portrait,
% 4. draws describing-function plots similar to Slotine Fig. 5.24 and 5.25.

clear; clc; close all;

%% Parameters for Numerical Simulation

% Van der Pol parameter
mu = 1; % It is well known that this system has a unique limit cycle for large mu

% Simulation time
tspan = [0 60]; % for a 60 seconds run

% Initial conditions:


x0 = [0.5 0];   % x(0) = 0.5 and dx/dt(0) = 0

%% Numerical simulation of the The Van der Pol Oscillator

% The Van der Pol equation is
%
%   xddot + mu*(x^2 - 1)*xdot + x = 0
%
% State variables:
% x(1) = x
% x(2) = dx/dt
[t,x] = ode45(@(t,x) [x(2); mu*(1-x(1)^2)*x(2)-x(1)],tspan,x0);

% Steady-state part: The response after 40 seconds is used to reduce the effect of the initial transient.
i = t >= 40;        % steady-state indexes of stamps
t_ss = t(i);        % steady state time stamps (crossing times: finding the time of each complete cycle)
x_ss = x(i,1);      % states at steady-state to be used for drawing the phase portraits. x(i,1) = x(1) = x in ode45 solver function
xdot_ss = x(i,2);   % derivatives of states at steady-state to be used for sketching the phase portraits.  x(i,2) = x(2) = xdot in ode45 solver function 

% Numerical amplitude: For a nearly symmetric oscillation: A = (maximum value - minimum value)/2
A_num = (max(x_ss)-min(x_ss))/2;

% to estimate the amplitude of the numerical limit cycle from the
% steady-state simulation for a symmetric periodic oscillation,
% (A - (-A)) / 2 = 2 A / 2 = A; (we can say that it's the mean of positive values and negative values)

% Numerical frequency: Find positive-going zero crossings of x(t).
z = find(x_ss(1:end-1)<0 & x_ss(2:end)>=0);

% The key point is:
% We do not need z to find the limit-cycle amplitude. We need it z to identify complete oscillation cycles and calculate the numerical period/frequency. 
% The reason for specifically choosing positive-going zero crossings instead of every zero crossing is that 
% using all zero crossings would give approximately half of the period, 
% because a sinusoidal-type oscillation crosses zero twice per cycle

% Consecutive zero crossings give approximately one period.
T_num = mean(diff(t_ss(z)));      % taking the average period

% Convert the period into frequency in Hz and rad/s.
f_num = 1/T_num;        % frequency in Hz
w_num = 2*pi*f_num;     % frequency in rad/s

% Display the numerical limit-cycle results.
fprintf('\nNumerical limit cycle:\n');
fprintf('A = %.4f\n', A_num);
fprintf('w = %.4f rad/s\n', w_num); % w is found from the numerical period that comes from the average of positive and negative crossings period
fprintf('f = %.4f Hz\n', f_num);
fprintf('T = %.4f s\n\n', T_num);

% Describing-function prediction for this example.
fprintf('Describing function prediction:\n');
fprintf('A = 2 \n');
fprintf('w = 1 rad/s\n');

%% Time response of Van Der Pol Oscillator
% Plot x(t) for the complete simulation.
figure(Name='Time Response', NumberTitle="off")
plot(t, x(:,1),'LineWidth', 1.5);     % plot states, x
grid on;
xlabel('Time (s)');
ylabel('x(t)');
title('Time Response','Van der Pol Oscillator');

%% Phase portrait of Van Der Pol Oscillator for the Limit Cycle Analysis
% The phase portrait shows dx/dt versus x.
figure(Name='Phase Portrait', NumberTitle="off")

plot(x(:,1), x(:,2),'LineWidth', 1);    % Inside: (A < 2), Transient-State Limit Cycle, Trajectory: Plotting x versus xdot in phase portrait 

hold on;

% Plot the steady-state limit cycle with a thicker line.
plot(x_ss, xdot_ss, 'LineWidth',2); % Outside: (A > 2), Steady-State Limit Cycle, ~ Numerical Limit Cycle

% The describing-function prediction is A = 2 and w = 1 rad/s.
% Therefore:
% x = 2*sin(theta)
% dx/dt = 2*cos(theta)

numpts = 360;     % number of points taken for resolution (here, 360 points are taken for each degrees)

theta = linspace(0, 2*pi, numpts);

plot(2*sin(theta), 2*cos(theta), 'k--', 'LineWidth', 1.5); % (A = 2), Theoretical Limit Cycle, DF Prediction 

grid on;
axis equal;
xlabel('x');
ylabel('dx/dt');
title('Phase Portrait','Van Der Pol Oscillator');
legend('Inside (A < 2): Transient-State Limit Cycle (Trajectory)', 'Outside (A > 2), Steady-State Limit Cycle, (~= Numerical Limit Cycle)', 'Onside (A = 2), Theoretical Limit Cycle, (DF Prediction)','Location','Best');

%% Describing function Generation
% The equations used in the describing-function analysis are:
% G(s) = 1/(s^2-mu*s+1)
% N(A,w) = j*mu*w*A^2/4
% Limit cycle: G(jw)N(A,w) = -1

% Frequency vector used for the complex-plane plots.
w = linspace(0.05,3,1000);

% Linear frequency response G(jw).
G = 1./((1-w.^2)-1i*mu*w);

%% Plot similar to Slotine's Book Fig. 5.24 Limjt cycle detection for frequency-dependent describing functions for different frequencies
% This plot compares G(jw) with a family of -1/N(A,w) curves.
figure(Name='Slotine''s Book Fig. 5.24', NumberTitle="off")

plot(real(G), imag(G), 'LineWidth', 2);

hold on;

% Amplitude range for the describing-function curves.
A = linspace(0.5,4,400);

% Fixed frequency values used to draw the family of -1/N curves.
w0 = [0.5 0.8 1 1.5 2];

for k = 1:length(w0)
    % Describing function for the current fixed frequency.
    N = 1i*mu*w0(k)*A.^2/4;

    % Negative inverse describing function.
    M = -1./N;

    % Plot in the complex plane.
    plot(real(M), imag(M), '--', 'LineWidth', 1, 'DisplayName', string(w0(k)));
end

% Mark the describing-function solution A = 2, w = 1 rad/s.
plot(0, 1/mu, 'ko', 'MarkerFaceColor', 'k');
text(0.05, 1/mu + 0.5, 'A = 2, \omega = 1');

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Nyquist Plot');
txt1 = {'Van der Pol Oscillator Limit Cycle Detection for different \omega, frequencies','(Slotine''s Book Fig. 5.24 Finding frequency-dependent DFs for different frequencies)'};
subtitle(txt1);
legend('G(j\omega), Linear TF', '-1/N(A, \omega) at \omega_{1} = 0.5', '-1/N(A, \omega) at \omega_{2} = 0.8', '-1/N(A, \omega) at \omega_{3} = 1', '-1/N(A, \omega) at \omega_{4} = 1.5', '-1/N(A, \omega) at \omega_{5} = 2', 'Limit Cycle');

%% Plot similar to Slotine's Book Fig. 5.25 Limjt cycle detection for frequency-dependent describing functions for different amplitudes
% This plot shows G(jw)N(A,w) for several fixed values of A.
% A limit cycle is predicted when a curve passes through (-1,0).
figure(Name='Slotine''s Book Fig. 5.25', NumberTitle="off")
hold on;

% Fixed amplitude values.
A_values = [1 1.5 2 2.5 3];

result = {};

for k = 1:length(A_values)
    % Select the current amplitude.
    A = A_values(k);

    % Calculate N(A,w).
    N = 1i*mu*w*A^2/4;

    % Calculate G(jw)N(A,w).
    L = G.*N;

    % Plot the result in the complex plane.
    plot(real(L), imag(L), 'LineWidth', 1.5);
        
end

% Mark the critical point (-1,0).
plot(-1,0, 'ko','MarkerFaceColor','k');
text(-0.95, 0.05, '(-1,0): Limit Cycle');

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Nyquist Plot'); 
txt2={'Van der Pol Oscillator Limit Cycle Detection for different amplitudes','(Slotine''s Book Fig. 5.24 Finding frequency-dependent describing functions for different amplitudes)'};
subtitle(txt2);
legend('A = 1','A = 1.5','A = 2','A = 2.5','A = 3','(-1,0): Limit Cycle Crossing');

