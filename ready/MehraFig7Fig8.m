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

figure;
hold on;
grid on;
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

        % Reconstruct theta
        theta = 6.02372 * x_out(:, 2) + 7.346 * x_out(:, 3);

        % Extract Steady-State Amplitude (evaluate the last 30% of the simulation)
        idx_ss = find(t_out > tspan(end) * 0.7);
        if ~isempty(idx_ss)
            theta_ss = theta(idx_ss);

            % Peak-to-Peak Amplitude / 2
            amplitude = (max(theta_ss) - min(theta_ss)) / 2;

            % Threshold to filter out numerical noise from stable equilibrium points
            if amplitude < 0.05
                amplitude = 0; 
            end

            amp_vec(j) = amplitude;
        end
    end

    % Plot Kp vs Amplitude for the current tau
    plot(Kp_vec, amp_vec, '-o', 'Color', colors(i), 'LineWidth', 1.5, ...
        'DisplayName', ['\tau = ' num2str(tau) ' s']);
end

% Formatting
xlabel('Pilot Gain (K_p)');
ylabel('Steady-State Limit Cycle Amplitude \theta (deg)');
title('Pilot Gain vs. Limit Cycle Amplitude (\theta) for Various Delays');
legend('Location', 'northwest');