function PIO_LimitCycleApp()
    % Global flag to handle Ctrl+C / Stop behavior
    stopRequested = false;

    % Create the main UI Figure and bind the close event to the stop function
    fig = uifigure('Name', 'PIO Limit Cycle & Phase Portrait Analysis (Mehra''s 1998 Paper)', ...
                   'Position', [100, 100, 1000, 800], ...
                   'CloseRequestFcn', @(src, event) stopSim());
    
    % Create a grid layout to separate inputs (left) from the plots (right)
    gl = uigridlayout(fig, [1 2]);
    gl.ColumnWidth = {340, '1x'}; 
    
    % --- Left Panel: Input Parameters ---
    inputPanel = uipanel(gl, 'Title', 'Simulation Parameters');
    
    inputGrid = uigridlayout(inputPanel, [15 2]);
    inputGrid.RowHeight = repmat({30}, 1, 15);
    inputGrid.RowHeight{7} = 'fit'; 
    inputGrid.RowHeight{8} = 80;    
    inputGrid.RowHeight{15} = 40;   
    inputGrid.ColumnWidth = {'1x', '1x'};
    
    % Helper functions for UI fields 
    function field = createNumField(parent, row, labelText, defaultVal)
        lbl = uilabel(parent, 'Interpreter', 'latex', 'Text', labelText, 'FontSize', 14);
        lbl.Layout.Row = row; lbl.Layout.Column = 1;
        field = uieditfield(parent, 'numeric', 'Value', defaultVal);
        field.Layout.Row = row; field.Layout.Column = 2;
    end

    function field = createTxtField(parent, row, labelText, defaultVal)
        lbl = uilabel(parent, 'Interpreter', 'latex', 'Text', labelText, 'FontSize', 14);
        lbl.Layout.Row = row; lbl.Layout.Column = 1;
        field = uieditfield(parent, 'text', 'Value', defaultVal);
        field.Layout.Row = row; field.Layout.Column = 2;
    end

    % Define parameter input fields using inline LaTeX
    lbl_defaults = uilabel(inputGrid, 'Text', '--- System Defaults ---', 'FontWeight', 'bold');
    lbl_defaults.Layout.Row = 1; 
    lbl_defaults.Layout.Column = [1 2];

    f_K      = createNumField(inputGrid, 2, 'Gain ($K$):', 20);
    f_S      = createNumField(inputGrid, 3, 'Rate Limit Max ($S$):', 15);
    f_R      = createNumField(inputGrid, 4, 'Rate Limit Min ($R$):', -15);
    f_thetac = createNumField(inputGrid, 5, 'Command ($\theta_c$):', 1);
    f_x0     = createTxtField(inputGrid, 6, 'Initial Cond ($x_0$):', '0, 0, 0, 0');
    
    % LaTeX description for Initial Conditions
    lbl_x0_desc = uilabel(inputGrid, 'Interpreter', 'latex', 'FontSize', 12, 'WordWrap', 'on', ...
        'Text', '$x_{1}$: Elevator Deflection, $x_{2}$, $x_{3}$, $x_{4}$: Short Period Dynamics of the Longitudinal Axis');
    lbl_x0_desc.Layout.Row = 7;
    lbl_x0_desc.Layout.Column = [1 2];
    
    % LaTeX Equation for the Saturation / Rate Limiter function
    eq_text = ['$$\dot{y} = \mathrm{SAT}(Ke) = \left\{ \begin{array}{ll} ' ...
               'S & \mathrm{if\ } Ke \ge S \\ ' ...
               'Ke & \mathrm{if\ } R < Ke < S \\ ' ...
               'R & \mathrm{if\ } Ke \le R \end{array} \right.$$'];
    lbl_eq = uilabel(inputGrid, 'Interpreter', 'latex', 'FontSize', 14, 'Text', eq_text);
    lbl_eq.Layout.Row = 8;
    lbl_eq.Layout.Column = [1 2];
    
    lbl_sweep = uilabel(inputGrid, 'Text', '--- Sweep Parameters ---', 'FontWeight', 'bold');
    lbl_sweep.Layout.Row = 9;
    lbl_sweep.Layout.Column = [1 2];

    f_Kp     = createTxtField(inputGrid, 10, 'Pilot Gains ($K_p$):', '1:0.5:15');
    f_tau    = createTxtField(inputGrid, 11, 'Delays ($\tau$):', '0, 0.03, 0.06, 0.09');
    f_tspan  = createTxtField(inputGrid, 12, 'Time Span ($t_{span}$):', '0, 60');
    
    lbl_target = uilabel(inputGrid, 'Text', '--- Phase Portrait Target ---', 'FontWeight', 'bold');
    lbl_target.Layout.Row = 13;
    lbl_target.Layout.Column = [1 2];
    
    % Phase portrait target inputs 
    targetGrid = uigridlayout(inputGrid, [1 4]);
    targetGrid.Layout.Row = 14; targetGrid.Layout.Column = [1 2];
    targetGrid.Padding = [0 0 0 0]; 
    targetGrid.ColumnWidth = {35, '1x', 30, '1x'};
    
    lbl_Kp_target = uilabel(targetGrid, 'Interpreter', 'latex', 'Text', '$K_{p}$:', 'FontSize', 14);
    lbl_Kp_target.Layout.Row = 1; lbl_Kp_target.Layout.Column = 1;
    f_Kp_target = uieditfield(targetGrid, 'numeric', 'Value', 10);
    f_Kp_target.Layout.Row = 1; f_Kp_target.Layout.Column = 2;
    
    lbl_tau_target = uilabel(targetGrid, 'Interpreter', 'latex', 'Text', '$\tau$:', 'FontSize', 14);
    lbl_tau_target.Layout.Row = 1; lbl_tau_target.Layout.Column = 3;
    f_tau_target = uieditfield(targetGrid, 'numeric', 'Value', 0.09);
    f_tau_target.Layout.Row = 1; f_tau_target.Layout.Column = 4;
    
    % --- Buttons ---
    % Run Button
    btnRun = uibutton(inputGrid, 'Text', 'Run Sweep', 'ButtonPushedFcn', @(btn,event) runSim());
    btnRun.Layout.Row = 15;
    btnRun.Layout.Column = 1;
    btnRun.BackgroundColor = [0 0.45 0.74];
    btnRun.FontColor = 'white';
    btnRun.FontWeight = 'bold';

    % Stop & Close Button
    btnStop = uibutton(inputGrid, 'Text', 'Stop & Close', 'ButtonPushedFcn', @(btn,event) stopSim());
    btnStop.Layout.Row = 15;
    btnStop.Layout.Column = 2;
    btnStop.BackgroundColor = [0.85 0.20 0.20]; 
    btnStop.FontColor = 'white';
    btnStop.FontWeight = 'bold';
    
    % --- Right Panel: Plotting Axes ---
    plotGrid = uigridlayout(gl, [2 1]);
    plotGrid.RowHeight = {'1x', '1x'};
    
    ax1 = uiaxes(plotGrid);
    title(ax1, 'Pilot Gain vs. Limit Cycle Amplitude (\theta) for Various Delays');
    xlabel(ax1, 'Pilot Gain (K_p)');
    ylabel(ax1, 'Limit Cycle Amplitude \theta (deg)');
    grid(ax1, 'on');
    
    ax2 = uiaxes(plotGrid);
    title(ax2, 'Steady-State Phase Portrait (\theta vs. d\theta/dt)');
    xlabel(ax2, '\theta (deg)');
    ylabel(ax2, 'Pitch Rate d\theta/dt (deg/s)');
    grid(ax2, 'on');
    
    % --- Automatic Start Timer ---
    % Wait 5 seconds, then automatically trigger the simulation
    t = timer('StartDelay', 1, 'ExecutionMode', 'singleShot', ...
              'TimerFcn', @(~,~) autoStart());
    start(t);
    
    function autoStart()
        if isvalid(fig)
            runSim();
        end
    end

    % --- Core Functions ---
    function stopSim()
        % Triggers the halt condition in the solvers and closes the app
        stopRequested = true;
        
        % Clean up the timer if the app is closed before 5 seconds
        timers = timerfindall;
        if ~isempty(timers)
            stop(timers); delete(timers);
        end
        
        if isvalid(fig)
            delete(fig);
        end
    end

    function status = checkStopFcn(~, ~, ~)
        % Output function injected into ode45/dde23 to allow Ctrl+C style interruption
        if stopRequested || ~isvalid(fig)
            status = 1; % Instructs the solver to halt immediately
        else
            status = 0;
        end
        drawnow limitrate; % Process UI button clicks during the loop
    end

    function runSim()
        stopRequested = false; % Reset stop flag on new run
        
        K = f_K.Value; S = f_S.Value; R = f_R.Value; thetac = f_thetac.Value;
        x0 = str2num(f_x0.Value); x0 = x0(:); %#ok<ST2NM>
        Kp_vec = str2num(f_Kp.Value); %#ok<ST2NM>
        tau_vec = str2num(f_tau.Value); %#ok<ST2NM>
        tspan = str2num(f_tspan.Value); %#ok<ST2NM>
        
        Kp_tgt = f_Kp_target.Value;
        tau_tgt = f_tau_target.Value;
        
        cla(ax1); hold(ax1, 'on');
        colors = lines(length(tau_vec));
        
        progDlg = uiprogressdlg(fig, 'Title', 'Running Sweep', 'Message', 'Integrating equations...');
        totalIters = length(tau_vec) * length(Kp_vec);
        currentIter = 0;
        
        % Attach the stop listener to the solvers
        opts_ode = odeset('OutputFcn', @checkStopFcn);
        opts_dde = ddeset('OutputFcn', @checkStopFcn);
        
        % 1. Execute Sweep for Top Plot
        for i = 1:length(tau_vec)
            if stopRequested || ~isvalid(fig), return; end
            
            tau = tau_vec(i);
            amp_vec = zeros(size(Kp_vec));
            
            for j = 1:length(Kp_vec)
                if stopRequested || ~isvalid(fig), return; end
                
                currentIter = currentIter + 1;
                if isvalid(progDlg), progDlg.Value = currentIter / totalIters; end
                
                Kp = Kp_vec(j);
                if tau == 0
                    sys_ode = @(t, x) [
                        max(R, min(S, K * (Kp * (thetac - (6.02372 * x(2) + 7.346 * x(3))) - x(1))));
                        x(3); x(4); x(1) - 5.29 * x(3) - 1.42 * x(4)];
                    [t_out, x_out] = ode45(sys_ode, tspan, x0, opts_ode);
                else
                    sys_dde = @(t, x, Z) [
                        max(R, min(S, K * (Kp * (thetac - (6.02372 * Z(2,1) + 7.346 * Z(3,1))) - x(1))));
                        x(3); x(4); x(1) - 5.29 * x(3) - 1.42 * x(4)];
                    sol = dde23(sys_dde, tau, x0, tspan, opts_dde);
                    if ~isempty(sol)
                        t_out = sol.x'; x_out = sol.y';
                    end
                end
                
                if stopRequested || ~isvalid(fig), return; end
                
                % Reconstruct theta and find amplitude
                theta = 6.02372 * x_out(:, 2) + 7.346 * x_out(:, 3);
                idx_ss = find(t_out > tspan(end) * 0.7);
                if ~isempty(idx_ss)
                    theta_ss = theta(idx_ss);
                    amplitude = (max(theta_ss) - min(theta_ss)) / 2;
                    if amplitude < 0.05, amplitude = 0; end
                    amp_vec(j) = amplitude;
                end
            end
            if isvalid(fig)
                plot(ax1, Kp_vec, amp_vec, '-o', 'Color', colors(i,:), ...
                    'LineWidth', 1.5, 'DisplayName', ['\tau = ' num2str(tau) ' s']);
            end
        end
        if isvalid(fig)
            legend(ax1, 'Location', 'best'); hold(ax1, 'off');
        end
        
        % 2. Execute Specific Run for Phase Portrait (Bottom Plot)
        if stopRequested || ~isvalid(fig), return; end
        if isvalid(progDlg), progDlg.Message = 'Generating Phase Portrait...'; end
        
        if tau_tgt == 0
            sys_ode_tgt = @(t, x) [
                max(R, min(S, K * (Kp_tgt * (thetac - (6.02372 * x(2) + 7.346 * x(3))) - x(1))));
                x(3); x(4); x(1) - 5.29 * x(3) - 1.42 * x(4)];
            [t_tgt, x_tgt] = ode45(sys_ode_tgt, tspan, x0, opts_ode);
        else
            sys_dde_tgt = @(t, x, Z) [
                max(R, min(S, K * (Kp_tgt * (thetac - (6.02372 * Z(2,1) + 7.346 * Z(3,1))) - x(1))));
                x(3); x(4); x(1) - 5.29 * x(3) - 1.42 * x(4)];
            sol_tgt = dde23(sys_dde_tgt, tau_tgt, x0, tspan, opts_dde);
            if ~isempty(sol_tgt)
                t_tgt = sol_tgt.x'; x_tgt = sol_tgt.y';
            end
        end
        
        if stopRequested || ~isvalid(fig), return; end
        
        theta_tgt = 6.02372 * x_tgt(:, 2) + 7.346 * x_tgt(:, 3);
        theta_dot_tgt = 6.02372 * x_tgt(:, 3) + 7.346 * x_tgt(:, 4);
        
        idx_phase = t_tgt > tspan(end) * 0.6;
        
        if isvalid(fig)
            cla(ax2);
            if any(idx_phase)
                plot(ax2, theta_tgt(idx_phase), theta_dot_tgt(idx_phase), 'k-', 'LineWidth', 1.5);
                title(ax2, sprintf('Steady-State Phase Portrait ($K_p$ = %g, $\\tau$ = %g s)', Kp_tgt, tau_tgt), 'Interpreter', 'latex', 'FontSize', 12);
            else
                title(ax2, 'Steady-State Phase Portrait (No Data)');
            end
        end
        
        if isvalid(progDlg)
            close(progDlg);
        end
    end
end