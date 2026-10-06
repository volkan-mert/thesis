% Bifurcation Analysis of Longitudinal Dynamics with Rate-Limited Actuator
clear; clc; close all;

% --- Fixed Parameters ---
K = 20;             % Actuator gain (from previous context)
theta_c = 0;        % Command input (zero for stability/bifurcation analysis)
sat_limit = 15;     % Saturation limit [-15, 15]

% --- Bifurcation Parameter Range ---
% Proportional Gain (Kp) sweep to detect onset of limit cycles (PIO)
Kp_min = 0.1;
Kp_max = 5.0;       % Adjust upper bound if necessary to see full bifurcation branch
Kp_steps = 300;
Kp_values = linspace(Kp_min, Kp_max, Kp_steps);

% --- Simulation Settings ---
tspan = [0 200];      % Total simulation time per Kp to ensure steady state
t_transient = 150;    % Discard early data to isolate steady-state/limit cycle behavior

% Initial conditions [x1, x2, x3, x4]
% x1: Actuator state (delta_e)
% x2, x3, x4: Airframe states
x0 = [0; 0; 0; 0]; 

% Prepare Figure
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;
xlabel('Proportional Gain ($K_p$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Steady-state Extrema of Actuator State ($x_1$)', 'Interpreter', 'latex', 'FontSize', 12);
title('Bifurcation Diagram: Actuator Rate Limiting', 'FontSize', 14);

% --- Main Bifurcation Loop ---
for i = 1:length(Kp_values)
    Kp = Kp_values(i);
    
    % Define ODE function handle
    ode_func = @(t, x) longitudinal_dynamics(t, x, Kp, K, sat_limit, theta_c);
    
    % Solve ODE using tight tolerances for accurate limit cycle detection
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [t, x] = ode45(ode_func, tspan, x0, options);
    
    % Isolate the steady-state portion
    steady_idx = t > t_transient;
    x1_steady = x(steady_idx, 1);
    
    if isempty(x1_steady)
        continue;
    end
    
    % Find local maxima to capture oscillatory behavior (limit cycles)
    [pks, ~] = findpeaks(x1_steady);
    
    if isempty(pks)
        % No peaks detected: System is at a stable equilibrium (fixed point)
        plot(Kp, x1_steady(end), 'b.', 'MarkerSize', 4);
    else
        % Peaks detected: System is in a limit cycle (oscillatory)
        plot(Kp * ones(size(pks)), pks, 'r.', 'MarkerSize', 4);
        
        % Optionally, you can also find and plot the valleys (minima) 
        % to show the full amplitude envelope of the limit cycle:
        [valleys, ~] = findpeaks(-x1_steady);
        plot(Kp * ones(size(valleys)), -valleys, 'r.', 'MarkerSize', 4);
    end
    
    % Update initial condition for the next step (continuation method)
    % This helps the solver track the solution branches smoothly
    x0 = x(end, :); 
end
hold off;

% ==========================================
%           STATE-SPACE DYNAMICS
% ==========================================
function dxdt = longitudinal_dynamics(~, x, Kp, K, sat_limit, theta_c)
    % Extract states
    x1 = x(1); % delta_e / y
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    
    % Output Equation
    % theta = 6.02372 * x2 + 7.346 * x3
    theta = 6.02372 * x2 + 7.346 * x3;
    
    % Rate-Limited Actuator Equation
    % v = K * (Kp * (theta_c - theta) - x1)
    v = K * (Kp * (theta_c - theta) - x1);
    
    % Saturation block [-15, 15]
    sat_v = max(-sat_limit, min(sat_limit, v)); 
    
    % --- State Derivatives ---
    dx1_dt = sat_v;
    dx2_dt = x3;
    dx3_dt = x4;
    dx4_dt = x1 - 5.29 * x3 - 1.42 * x4;
    
    % Return column vector
    dxdt = [dx1_dt; dx2_dt; dx3_dt; dx4_dt];
end