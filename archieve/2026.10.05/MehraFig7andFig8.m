% --- Pilot Gain vs. Limit Cycle Amplitude ---
% Parameters
K = 20;
S = 15;
R = -15;
thetac = 1; % 1-degree step command

% Initial Conditions [x1; x2; x3; x4]
x0 = [0; 0; 0; 0];
tspan = [0 60]; % Increased time to ensure steady-state limit cycle is reached

% Arrays for Analysis
tau_vec = [0, 0.03, 0.06, 0.09];
Kp_vec = 1:0.5:15; % Range of Pilot Gains to test

% Setup Main Figure
figure('Name', 'Pilot Model Analysis', 'Position', [100, 100, 1200, 500]);

% Setup Subplot 1 (Amplitude Plot)
ax1 = subplot(1,2,1);
hold(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'Pilot Gain ($K_p$)', 'Interpreter', 'latex');
ylabel(ax1, 'Steady-State Limit Cycle Amplitude $\theta$ (deg)', 'Interpreter', 'latex');
title(ax1, 'Figure 7 of Mehra 1998 Bifurcation and LC Analysis Paper', 'Interpreter', 'latex');

% Setup Subplot 2 (Phase Portrait)
ax2 = subplot(1,2,2);
hold(ax2, 'on');
grid(ax2, 'on');
xlabel(ax2, '$\theta$ (deg)', 'Interpreter', 'latex');
ylabel(ax2, '$\dot{\theta}$ (deg/s)', 'Interpreter', 'latex');
title(ax2, '$\theta$ vs. $\dot{\theta}$ Figure 8 of Mehra 1998 Bifurcation and LC Analysis Paper', 'Interpreter', 'latex');

colors = ['b', 'r', 'g', 'm'];

for i = 1:length(tau_vec)
    tau = tau_vec(i);
    amp_vec = zeros(size(Kp_vec));
    
    for j = 1:length(Kp_vec)
        Kp = Kp_vec(j);
        
        if tau == 0
            % Standard ODE without delay (dde23 requires lag > 0)
            sys_ode = @(t, x) [
                max(R, min(S, K * (Kp * (thetac - (6.02372 * x(2) + 7.346 * x(3))) - x(1))));
                x(3);
                x(4);
                x(1) - 5.29 * x(3) - 1.42 * x(4)
                ];
            [t_out, x_out] = ode45(sys_ode, tspan, x0);
        else
            % DDE with delay for pilot model Kp * exp(-s*tau)
            % Z(:,1) represents the state vector at t - tau
            sys_dde = @(t, x, Z) [
                max(R, min(S, K * (Kp * (thetac - (6.02372 * Z(2,1) + 7.346 * Z(3,1))) - x(1))));
                x(3);
                x(4);
                x(1) - 5.29 * x(3) - 1.42 * x(4)
                ];
            % Solve DDE using constant history vector x0
            sol = dde23(sys_dde, tau, x0, tspan);
            t_out = sol.x';
            x_out = sol.y';
        end
        
        % Reconstruct theta and theta_dot
        % Since x(2)' = x(3) and x(3)' = x(4)
        theta = 6.02372 * x_out(:, 2) + 7.346 * x_out(:, 3);
        theta_dot = 6.02372 * x_out(:, 3) + 7.346 * x_out(:, 4);
        
        % Extract Steady-State Amplitude (evaluate the last 30% of the simulation)
        idx_ss = find(t_out > tspan(end) * 0.7);
        
        if ~isempty(idx_ss)
            theta_ss = theta(idx_ss);
            theta_dot_ss = theta_dot(idx_ss);
            
            % Peak-to-Peak Amplitude / 2
            amplitude = (max(theta_ss) - min(theta_ss)) / 2;
            
            % Threshold to filter out numerical noise from stable equilibrium points
            if amplitude < 0.05
                amplitude = 0; 
            end
            amp_vec(j) = amplitude;
            
            % Plot Phase Portrait on Subplot 2
            % Only add the largest limit cycle to the legend to avoid clutter
            if j == length(Kp_vec) 
                plot(ax2, theta_ss, theta_dot_ss, 'Color', colors(i), 'DisplayName', ['$\tau = ' num2str(tau) '$ s']);
            elseif amplitude > 0 % Only draw the actual limit cycles, ignoring stable points
                plot(ax2, theta_ss, theta_dot_ss, 'Color', colors(i), 'HandleVisibility', 'off');
            end
        end
    end
    
    % Plot Kp vs Amplitude for the current tau on Subplot 1
    plot(ax1, Kp_vec, amp_vec, '-o', 'Color', colors(i), 'LineWidth', 1.5, ...
        'DisplayName', ['$\tau = ' num2str(tau) '$ s']);
end

% Formatting Legends
legend(ax1, 'Location', 'northwest', 'Interpreter', 'latex');
legend(ax2, 'Location', 'best', 'Interpreter', 'latex');