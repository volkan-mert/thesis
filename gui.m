function aircraft_limit_cycle_gui()
%% =========================================================================
%  INTERACTIVE LONGITUDINAL DYNAMICS & LIMIT CYCLE SIMULATOR
%  Includes Real-Time Simulation, DDE Sweep, and Live Nichols Chart
% =========================================================================

    % --- 1. SYSTEM CONSTANTS & AIRCRAFT DYNAMICS ---
    K = 20.0;           % Actuator gain
    
    % --- 2. INTERACTIVE SIMULATION PARAMETERS ---
    K_p   = 5.0;        % Pilot gain
    R_max = 15.0;       % Actuator saturation limit S & R [deg/s]
    tau   = 0.0;        % Pilot delay [s]
    thetac= 1.0;        % Step command [deg]
    
    % Simulation state flags
    is_running = true;
    sim_time   = 0.0;
    dt         = 0.002;  % Integration time step (500 Hz)
    
    % Aircraft state vector: [x1; x2; x3; x4]
    x1 = 0.0; x2 = 0.0; x3 = 0.0; x4 = 0.0;
    
    % Delay & Data history buffers 
    hist_len   = 500;
    theta_hist = zeros(1, hist_len);
    hist_idx   = 1;

    buf_len      = 7500;
    t_buf        = nan(1, buf_len);
    theta_buf    = nan(1, buf_len);
    cmd_buf      = nan(1, buf_len);
    tdot_buf     = nan(1, buf_len); 
    rdot_buf     = nan(1, buf_len);
    e_act_buf    = nan(1, buf_len); % Buffer for describing function input amplitude (Ai)

    % Nichols Chart Caching States (prevents lag by only redrawing when Ai changes)
    last_Ai = -1; last_Kp = -1; last_Rmax = -1;

    % --- 3. GUI FIGURE & CONTROL PANEL CREATION ---
    fig = figure('Name', 'Pilot-in-the-Loop Limit Cycle Simulator', ...
                 'Color', [0.95 0.95 0.96], 'Position', [50, 50, 1400, 850], ...
                 'NumberTitle', 'off', 'CloseRequestFcn', @close_callback);

    % TOP ROW: Real-Time Simulation Plots
    ax_phase = axes('Parent', fig, 'Position', [0.05, 0.62, 0.27, 0.32]);
    grid(ax_phase, 'on'); hold(ax_phase, 'on');
    h_phase_line = plot(ax_phase, 0, 0, 'Color', [0.85 0.32 0.10], 'LineWidth', 1.6);
    h_phase_curr = plot(ax_phase, 0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    xlabel(ax_phase, 'Pitch Angle \theta (deg)', 'FontWeight', 'bold');
    ylabel(ax_phase, 'Pitch Rate d\theta/dt (deg/s)', 'FontWeight', 'bold');
    title(ax_phase, 'Real-Time Phase Portrait', 'FontSize', 11);

    ax_theta = axes('Parent', fig, 'Position', [0.37, 0.62, 0.27, 0.32]);
    grid(ax_theta, 'on'); hold(ax_theta, 'on');
    h_cmd_line   = plot(ax_theta, 0, 0, 'r--', 'LineWidth', 1.2);
    h_theta_line = plot(ax_theta, 0, 0, 'Color', [0.00 0.45 0.74], 'LineWidth', 1.4);
    ylabel(ax_theta, '\theta (deg)', 'FontWeight', 'bold');
    title(ax_theta, 'Pitch Angle Response & Command', 'FontSize', 11);

    ax_rdot = axes('Parent', fig, 'Position', [0.69, 0.62, 0.27, 0.32]);
    grid(ax_rdot, 'on'); hold(ax_rdot, 'on');
    h_rdot_line = plot(ax_rdot, 0, 0, 'Color', [0.47 0.67 0.19], 'LineWidth', 1.4);
    h_lim_pos   = yline(ax_rdot, R_max, 'r--', 'LineWidth', 1.4);
    h_lim_neg   = yline(ax_rdot, -R_max, 'r--', 'LineWidth', 1.4);
    xlabel(ax_rdot, 'Simulation Time (s)', 'FontWeight', 'bold');
    ylabel(ax_rdot, 'dy/dt (deg/s)', 'FontWeight', 'bold');
    title(ax_rdot, 'Actuator Rate Saturation', 'FontSize', 11);

    % MIDDLE ROW: DDE Bifurcation Sweep & Nichols Plots
    ax_amp = axes('Parent', fig, 'Position', [0.05, 0.24, 0.27, 0.28]);
    grid(ax_amp, 'on'); hold(ax_amp, 'on');
    xlabel(ax_amp, 'Pilot Gain (K_p)', 'FontWeight', 'bold');
    ylabel(ax_amp, 'Limit Cycle Amplitude (deg)', 'FontWeight', 'bold');
    title(ax_amp, 'Pilot Gain vs. Steady-State Amplitude', 'FontSize', 11);
    text(ax_amp, 0.5, 0.5, 'Click "Compute Sweep" to generate', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'Color', [0.5 0.5 0.5]);

    ax_freq = axes('Parent', fig, 'Position', [0.37, 0.24, 0.27, 0.28]);
    grid(ax_freq, 'on'); hold(ax_freq, 'on');
    xlabel(ax_freq, 'Pilot Gain (K_p)', 'FontWeight', 'bold');
    ylabel(ax_freq, 'Limit Cycle Frequency (Hz)', 'FontWeight', 'bold');
    title(ax_freq, 'Pilot Gain vs. Limit Cycle Frequency', 'FontSize', 11);
    text(ax_freq, 0.5, 0.5, 'Click "Compute Sweep" to generate', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'Color', [0.5 0.5 0.5]);

    % Target Axes for Nichols Chart
    ax_nichols = axes('Parent', fig, 'Position', [0.69, 0.24, 0.27, 0.28]);

    % --- 4. INTERACTIVE CONTROLS (PANEL AT BOTTOM) ---
    panel = uipanel('Parent', fig, 'Position', [0.02, 0.02, 0.96, 0.14], 'BackgroundColor', 'w');

    uicontrol('Parent', panel, 'Style', 'text', 'Position', [20, 60, 200, 20], 'String', 'Actuator Saturation S (deg/s):', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'BackgroundColor', 'w');
    lbl_rmax = uicontrol('Parent', panel, 'Style', 'text', 'Position', [230, 60, 60, 20], 'String', sprintf('%.1f', R_max), 'HorizontalAlignment', 'left', 'ForegroundColor', 'b', 'BackgroundColor', 'w');
    sld_rmax = uicontrol('Parent', panel, 'Style', 'slider', 'Position', [20, 40, 250, 22], 'Min', 5, 'Max', 30, 'Value', R_max, 'Callback', @update_rmax);

    uicontrol('Parent', panel, 'Style', 'text', 'Position', [320, 60, 200, 20], 'String', 'Pilot Gain K_p:', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'BackgroundColor', 'w');
    lbl_kp = uicontrol('Parent', panel, 'Style', 'text', 'Position', [530, 60, 60, 20], 'String', sprintf('%.2f', K_p), 'HorizontalAlignment', 'left', 'ForegroundColor', 'b', 'BackgroundColor', 'w');
    sld_kp = uicontrol('Parent', panel, 'Style', 'slider', 'Position', [320, 40, 250, 22], 'Min', 0.5, 'Max', 15.0, 'Value', K_p, 'Callback', @update_kp);

    uicontrol('Parent', panel, 'Style', 'text', 'Position', [620, 60, 200, 20], 'String', 'Pilot Delay \tau (s):', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'BackgroundColor', 'w');
    lbl_tau = uicontrol('Parent', panel, 'Style', 'text', 'Position', [830, 60, 60, 20], 'String', sprintf('%.2f', tau), 'HorizontalAlignment', 'left', 'ForegroundColor', 'b', 'BackgroundColor', 'w');
    sld_tau = uicontrol('Parent', panel, 'Style', 'slider', 'Position', [620, 40, 250, 22], 'Min', 0.0, 'Max', 0.2, 'Value', tau, 'Callback', @update_tau);

    btn_pause   = uicontrol('Parent', panel, 'Style', 'pushbutton', 'Position', [20, 5, 100, 30], 'String', 'Pause', 'FontWeight', 'bold', 'Callback', @toggle_pause);
    btn_perturb = uicontrol('Parent', panel, 'Style', 'pushbutton', 'Position', [130, 5, 180, 30], 'String', 'Toggle Command (\theta_c)', 'FontWeight', 'bold', 'Callback', @toggle_command);
    btn_reset   = uicontrol('Parent', panel, 'Style', 'pushbutton', 'Position', [320, 5, 100, 30], 'String', 'Reset State', 'FontWeight', 'bold', 'Callback', @reset_sim);
    btn_sweep   = uicontrol('Parent', panel, 'Style', 'pushbutton', 'Position', [430, 5, 220, 30], 'String', 'Compute Bifurcation Sweep', 'FontWeight', 'bold', 'BackgroundColor', [0.8 1 0.8], 'Callback', @run_sweep);
    btn_stop    = uicontrol('Parent', panel, 'Style', 'pushbutton', 'Position', [660, 5, 120, 30], 'String', 'Stop & Close', 'FontWeight', 'bold', 'BackgroundColor', [1 0.6 0.6], 'Callback', @close_callback);
    status_box  = uicontrol('Parent', panel, 'Style', 'text', 'Position', [790, 5, 200, 30], 'String', 'Status: Running...', 'FontWeight', 'bold', 'ForegroundColor', [0 0.5 0], 'BackgroundColor', 'w', 'FontSize', 11);

    % --- 5. REAL-TIME SIMULATION & ANIMATION LOOP ---
    step_skip = 12; % Integrate 12 steps before drawing to maintain UI responsiveness
    t_window  = 15; % Seconds of history shown on the graphs
    
    while ishandle(fig)
        if is_running
            for s = 1:step_skip
                theta = 6.02372 * x2 + 7.346 * x3;
                theta_hist(hist_idx) = theta;
                
                % Extract delayed theta
                delay_steps = round(tau / dt);
                del_idx = hist_idx - delay_steps;
                if del_idx < 1
                    del_idx = del_idx + hist_len;
                end
                theta_del = theta_hist(del_idx);
                
                % Controller law & actuator error
                u = K_p * (thetac - theta_del);
                e_act = u - x1;
                
                % Actuator state derivative with saturation
                x1_dot_raw = K * e_act;
                if x1_dot_raw > R_max
                    x1_dot = R_max;
                elseif x1_dot_raw < -R_max
                    x1_dot = -R_max;
                else
                    x1_dot = x1_dot_raw;
                end
                
                % Aircraft Phase-Variable Dynamics
                x2_dot = x3;
                x3_dot = x4;
                x4_dot = x1 - 5.29 * x3 - 1.42 * x4;
                
                % Numerical integration (Euler forward)
                x1 = x1 + dt * x1_dot; x2 = x2 + dt * x2_dot;
                x3 = x3 + dt * x3_dot; x4 = x4 + dt * x4_dot;
                
                sim_time = sim_time + dt;
                
                hist_idx = hist_idx + 1;
                if hist_idx > hist_len, hist_idx = 1; end
            end
            
            theta     = 6.02372 * x2 + 7.346 * x3;
            theta_dot = 6.02372 * x3 + 7.346 * x4;
            
            % Push to circular display buffer
            t_buf        = [t_buf(2:end), sim_time];
            theta_buf    = [theta_buf(2:end), theta];
            cmd_buf      = [cmd_buf(2:end), thetac];
            tdot_buf     = [tdot_buf(2:end), theta_dot];
            rdot_buf     = [rdot_buf(2:end), x1_dot];
            e_act_buf    = [e_act_buf(2:end), e_act];
            
            % Update visual graphics
            set(h_phase_line, 'XData', theta_buf, 'YData', tdot_buf);
            set(h_phase_curr, 'XData', theta, 'YData', theta_dot);
            set(h_cmd_line, 'XData', t_buf, 'YData', cmd_buf);
            set(h_theta_line, 'XData', t_buf, 'YData', theta_buf);
            set(h_rdot_line, 'XData', t_buf, 'YData', rdot_buf);
            
            % Adaptive Scaling for Phase Portrait, Pitch Angle & Actuator Rate
            valid_theta = theta_buf(~isnan(theta_buf));
            valid_tdot  = tdot_buf(~isnan(tdot_buf));
            valid_rdot  = rdot_buf(~isnan(rdot_buf));
            
            if ~isempty(valid_theta) && ~isempty(valid_tdot)
                theta_bound = max(10, max(abs(valid_theta)) * 1.15);
                tdot_bound  = max(20, max(abs(valid_tdot)) * 1.15);
                xlim(ax_phase, [-theta_bound, theta_bound]);
                ylim(ax_phase, [-tdot_bound, tdot_bound]);
                ylim(ax_theta, [-theta_bound, theta_bound]);
            end
            
            if ~isempty(valid_rdot)
                rdot_bound = max(10, max(R_max, max(abs(valid_rdot))) * 1.25);
                ylim(ax_rdot, [-rdot_bound, rdot_bound]);
            end
            
            xlim(ax_theta, [max(0, sim_time - t_window), max(t_window, sim_time)]);
            xlim(ax_rdot, [max(0, sim_time - t_window), max(t_window, sim_time)]);

            % -------------------------------------------------------------
            % Live Nichols Chart with User Script Embedded
            % -------------------------------------------------------------
            valid_e_act = e_act_buf(~isnan(e_act_buf));
            if ~isempty(valid_e_act)
                Ai_val = max(abs(valid_e_act));
            else
                Ai_val = K_p; 
            end
            if Ai_val < 0.1, Ai_val = 0.1; end 
            
            % Update Nichols Chart only if parameters shift significantly (prevents GUI lag)
            if abs(Ai_val - last_Ai) > 0.05 || K_p ~= last_Kp || R_max ~= last_Rmax
                
                %% 1. Linear system G(s)
                Kp_df = K_p;       % Pilot gain (from GUI)
                M_del_e = 0.537;   % Input signal entering the rate limiter
                omega_n = 2.3;     % Natural frequency
                zeta_sp = 1.42 / omega_n / 2;   % Damping ratio
                num = Kp_df*M_del_e*[1 0.82];
                den = [1 1.42 omega_n^2 0];
                Gs = tf(num,den);
                
                %% 2. Frequency range for G(jw)
                w_G = linspace(0.3,3.5,1000);
                
                %% 3. Rate Limiter Element parameters
                R_df  = R_max;     % Rate limit, deg/s (from GUI)
                Ai = Ai_val;       % Input amplitude, deg (Dynamic from buffer)
                w_N = linspace(0.01, 10, 1001);   % Frequencies, rad/s
                w_onset = R_df / Ai;
                
                %% 4. Describing function N(Ai,w)
                N = zeros(size(w_N));
                for k = 1:length(w_N)
                    alpha = w_N(k) / w_onset;
                    % Region I
                    if alpha <= 1
                        M = 1;
                        phi = 0;
                    % Region II
                    elseif alpha < 1.862
                        M = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
                        phi = 0.528*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171;
                    % Region III
                    else
                        varpi = w_onset / w_N(k);
                        M = (4/pi)*varpi;
                        phi = -acos((pi/2)*varpi);
                    end
                    N(k) = M*exp(1j*phi);
                end
                
                %% 5. Calculate -1/N(Ai,w)
                minus_inv_N = -1 ./ N;
                
                %% 6. Convert -1/N to an FRD model
                response = reshape(minus_inv_N,1,1,[]);
                sys_N = frd(response,w_N);
                
                %% 7. Nichols Chart
                cla(ax_nichols);
                
                % Force colors to match the user's screenshot layout
                ax_nichols.ColorOrder = [0 0 1; 1 0 0]; % Blue for G, Red for sys_N
                
                p1 = nicholsplot(ax_nichols, Gs, w_G);
                hold(ax_nichols, 'on');
                p2 = nicholsplot(ax_nichols, sys_N);
                
                %% 8. Shift -1/N from +180 deg to -180 deg
                p2.PhaseMatchingEnabled = 'on';
                p2.PhaseMatchingFrequency = w_N(1);
                phase_first = rad2deg(angle(minus_inv_N(1)));
                p2.PhaseMatchingValue = phase_first - 360;
                
                %% 9. Figure settings
                grid(ax_nichols, 'on');
                % Limits removed to allow adaptive boundary framing as requested
                yline(ax_nichols, 0,'--k','DisplayName','0 dB');
                title(ax_nichols, 'Nichols Chart of G(j\omega) and -1/N(A_i,\omega)', 'FontSize', 11);
                legend(ax_nichols, 'G(j\omega)','-1/N(A_i,\omega)','0 dB', 'Location', 'northeast');
                
                % Cache state
                last_Ai = Ai_val;
                last_Kp = K_p;
                last_Rmax = R_max;
            end
        end
        drawnow limitrate;
    end

    % --- 6. NESTED CALLBACK FUNCTIONS ---
    function update_rmax(~, ~)
        R_max = get(sld_rmax, 'Value');
        set(lbl_rmax, 'String', sprintf('%.1f', R_max));
        set(h_lim_pos, 'Value', R_max); set(h_lim_neg, 'Value', -R_max);
    end

    function update_kp(~, ~)
        K_p = get(sld_kp, 'Value');
        set(lbl_kp, 'String', sprintf('%.2f', K_p));
    end

    function update_tau(~, ~)
        tau = get(sld_tau, 'Value');
        set(lbl_tau, 'String', sprintf('%.2f', tau));
    end

    function toggle_pause(~, ~)
        is_running = ~is_running;
        if is_running
            set(btn_pause, 'String', 'Pause');
            set(status_box, 'String', 'Status: Running...', 'ForegroundColor', [0 0.5 0]);
        else
            set(btn_pause, 'String', 'Resume');
            set(status_box, 'String', 'Status: Paused', 'ForegroundColor', [0.8 0 0]);
        end
    end

    function toggle_command(~, ~)
        thetac = -thetac; 
    end

    function reset_sim(~, ~)
        x1 = 0.0; x2 = 0.0; x3 = 0.0; x4 = 0.0;
        thetac = 1.0; sim_time = 0.0; hist_idx = 1;
        theta_hist(:) = 0.0; t_buf(:) = nan; theta_buf(:) = nan;
        cmd_buf(:) = nan; tdot_buf(:) = nan; rdot_buf(:) = nan; e_act_buf(:) = nan;
    end

    function run_sweep(~, ~)
        was_running = is_running;
        is_running = false; 
        
        h_wait = waitbar(0, 'Initializing Sweep...', 'Name', 'Limit Cycle Analysis');
        
        tau_vec = [0, 0.03, 0.06, 0.09];
        Kp_vec  = 1:0.5:15;
        tspan_sw = [0 60];
        x0_sw    = [0; 0; 0; 0];
        thetac_sw = 1; 
        
        cla(ax_amp); cla(ax_freq);
        colors = ['b', 'r', 'g', 'm'];
        total_iters = length(tau_vec) * length(Kp_vec);
        iter = 0;
        
        for i_tau = 1:length(tau_vec)
            tau_sw = tau_vec(i_tau);
            amp_vec = zeros(size(Kp_vec));
            freq_vec = zeros(size(Kp_vec));
            
            for j_kp = 1:length(Kp_vec)
                if ~ishandle(fig)
                    if ishandle(h_wait), close(h_wait); end
                    return; 
                end
                
                Kp_sw = Kp_vec(j_kp);
                iter = iter + 1;
                if ishandle(h_wait)
                    waitbar(iter/total_iters, h_wait, sprintf('Simulating \\tau = %.2f, K_p = %.1f (%d/%d)', tau_sw, Kp_sw, iter, total_iters));
                end
                
                if tau_sw == 0
                    sys_ode = @(t, x) [
                        max(R_max * -1, min(R_max, K * (Kp_sw * (thetac_sw - (6.02372 * x(2) + 7.346 * x(3))) - x(1))));
                        x(3); x(4);
                        x(1) - 5.29 * x(3) - 1.42 * x(4)
                    ];
                    [t_out, x_out] = ode45(sys_ode, tspan_sw, x0_sw);
                else
                    sys_dde = @(t, x, Z) [
                        max(R_max * -1, min(R_max, K * (Kp_sw * (thetac_sw - (6.02372 * Z(2,1) + 7.346 * Z(3,1))) - x(1))));
                        x(3); x(4);
                        x(1) - 5.29 * x(3) - 1.42 * x(4)
                    ];
                    sol = dde23(sys_dde, tau_sw, x0_sw, tspan_sw);
                    t_out = sol.x'; x_out = sol.y';
                end
                
                theta_out = 6.02372 * x_out(:, 2) + 7.346 * x_out(:, 3);
                idx_ss = find(t_out > tspan_sw(end) * 0.7);
                
                if ~isempty(idx_ss) && length(idx_ss) > 10
                    theta_ss = theta_out(idx_ss);
                    t_ss = t_out(idx_ss);
                    
                    amplitude = (max(theta_ss) - min(theta_ss)) / 2;
                    if amplitude < 0.05
                        amp_vec(j_kp) = 0; freq_vec(j_kp) = 0;
                    else
                        amp_vec(j_kp) = amplitude;
                        theta_cen = theta_ss - mean(theta_ss);
                        zc = find(theta_cen(1:end-1) .* theta_cen(2:end) < 0);
                        if length(zc) >= 2
                            period = 2 * mean(diff(t_ss(zc)));
                            freq_vec(j_kp) = 1 / period;
                        else
                            freq_vec(j_kp) = 0;
                        end
                    end
                end
            end
            plot(ax_amp, Kp_vec, amp_vec, '-o', 'Color', colors(i_tau), 'LineWidth', 1.5, 'DisplayName', ['\tau = ' num2str(tau_sw) ' s']);
            plot(ax_freq, Kp_vec, freq_vec, '-s', 'Color', colors(i_tau), 'LineWidth', 1.5, 'DisplayName', ['\tau = ' num2str(tau_sw) ' s']);
        end
        
        if ishandle(h_wait)
            close(h_wait);
        end
        legend(ax_amp, 'Location', 'northwest');
        legend(ax_freq, 'Location', 'northeast');
        is_running = was_running; 
    end

    function close_callback(~, ~)
        is_running = false;
        if ishandle(fig)
            delete(fig);
        end
    end
end