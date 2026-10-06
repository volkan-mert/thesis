% Bifurcation Analysis of Actuator-Plant System
clear; clc; close all;

% --- Fixed Parameters from Image ---
K = 20;
S = 15;
theta_c = 0; % Command input (assume 0 for stability/bifurcation analysis)

% --- Bifurcation Parameter Range ---
% Assuming Proportional Gain (Kp) is the bifurcation parameter
Kp_min = 0.1;
Kp_max = 10;
Kp_steps = 200;
Kp_values = linspace(Kp_min, Kp_max, Kp_steps);

% --- Simulation Settings ---
tspan = [0 150];      % Total simulation time per Kp
t_transient = 100;    % Time to discard to allow transients to decay

% Initial conditions [x1, x2, x3]
x0 = [0; 0; 0]; 

% Prepare Figure
figure('Color', 'w');
hold on; grid on;
xlabel('Proportional Gain (K_p)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Steady-state / Extrema of $x_1$', 'Interpreter', 'latex', 'FontSize', 12);
title('Bifurcation Diagram', 'FontSize', 14);

% --- Main Bifurcation Loop ---
for i = 1:length(Kp_values)
    Kp = Kp_values(i);
    
    % Define ODE 
    ode_func = @(t, x) system_dynamics(t, x, Kp, K, S, theta_c);
    
    % Solve ODE
    % Using somewhat tight tolerances to accurately capture limit cycles
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [t, x] = ode45(ode_func, tspan, x0, options);
    
    % Isolate the steady-state portion (discard transient)
    steady_idx = t > t_transient;
    x1_steady = x(steady_idx, 1);
    
    if isempty(x1_steady)
        continue;
    end
    
    % Find local maxima (peaks) to capture oscillatory behavior (limit cycles)
    [pks, ~] = findpeaks(x1_steady);
    
    if isempty(pks)
        % If no peaks, the system settled at a stable equilibrium (fixed point)
        plot(Kp, x1_steady(end), 'k.', 'MarkerSize', 5);
        
        % Update initial condition for the next step (continuation method)
        x0 = x(end, :); 
    else
        % If peaks exist, plot all unique peak values (indicates limit cycles/chaos)
        plot(Kp * ones(size(pks)), pks, 'k.', 'MarkerSize', 5);
        
        % Update initial condition
        x0 = x(end, :);
    end
end
hold off;

% ==========================================
%           SYSTEM DYNAMICS FUNCTION
% ==========================================
function dxdt = system_dynamics(t, x, Kp, K, S, theta_c)
    % State variables
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    
    % Actuator Equation (from provided image)
    % e = Kp * [theta_c - (6.02372*x2 + 7.346*x3)] - x1
    v = K * (Kp * (theta_c - 6.02372 * x2 - 7.346 * x3) - x1);
    
    % Saturation function SAT(v)
    sat_v = max(-S, min(S, v)); 
    
    % State 1 Derivative (Actuator)
    dx1_dt = sat_v;
    
    % ----------------------------------------------------
    % MISSING PLANT DYNAMICS: 
    % The image does not provide the equations for x2 and x3.
    % Replace the placeholder equations below with your actual 
    % aircraft dynamic model (e.g., short-period approximation).
    % ----------------------------------------------------
    
    % PLACEHOLDERS (Replace these lines)
    dx2_dt = -0.5 * x2 + x1; 
    dx3_dt = -0.1 * x3 + x2; 
    
    % Return column vector
    dxdt = [dx1_dt; dx2_dt; dx3_dt];
end