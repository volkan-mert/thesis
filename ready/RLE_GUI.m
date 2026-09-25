function RLE_GUI()
    % 1. Create the main UI figure window
    fig = uifigure('Name', 'Actuator Dynamics & Rate Limiting', 'Position', [100, 100, 950, 650]);

    % 2. Define the timer FIRST so it exists in the workspace for all callbacks
    autoRunTimer = timer('ExecutionMode', 'singleShot', ...
                         'StartDelay', 1, ...
                         'TimerFcn', @(~,~) runSim());

    % 3. Bind Window closing events (Standard 'X' and CTRL+C)
    fig.CloseRequestFcn = @(src, event) closeApp(src, autoRunTimer);
    fig.WindowKeyPressFcn = @(src, event) handleKeyPress(src, event, autoRunTimer);

    % Set up the main grid layout (left for controls, right for plots)
    gl = uigridlayout(fig, [2, 2]);
    gl.ColumnWidth = {280, '1x'};
    gl.RowHeight = {'1x', '1x'};

    % Control Panel Container
    pnl = uipanel(gl, 'Title', 'Actuator Parameters & Pilot Command');
    pnl.Layout.Row = [1 2];
    pnl.Layout.Column = 1;

    % Grid layout inside the control panel (Expanded to 9 rows for the equation)
    pGrid = uigridlayout(pnl, [9, 2]);
    pGrid.RowHeight = {30, 30, 30, 30, 30, 30, 40, 80, '1x'};
    pGrid.ColumnWidth = {'1x', '1x'};

    % --- UI Controls ---
    
    % Forward Gain (K)
    uilabel(pGrid, 'Text', 'Forward Gain (K):', 'HorizontalAlignment', 'right');
    editK = uieditfield(pGrid, 'numeric', 'Value', 20);

    % Upper Saturation (S)
    uilabel(pGrid, 'Text', 'Upper Rate Limit (S):', 'HorizontalAlignment', 'right');
    editS = uieditfield(pGrid, 'numeric', 'Value', 15);

    % Lower Saturation (R)
    uilabel(pGrid, 'Text', 'Lower Rate Limit (R):', 'HorizontalAlignment', 'right');
    editR = uieditfield(pGrid, 'numeric', 'Value', -15);

    % Input Amplitude
    uilabel(pGrid, 'Text', 'Command Amplitude (deg):', 'HorizontalAlignment', 'right');
    editAmp = uieditfield(pGrid, 'numeric', 'Value', 2);

    % Input Frequency
    uilabel(pGrid, 'Text', 'Command Frequency (rad/s):', 'HorizontalAlignment', 'right');
    editFreq = uieditfield(pGrid, 'numeric', 'Value', pi);

    % Simulation Time
    uilabel(pGrid, 'Text', 'Simulation Time (s):', 'HorizontalAlignment', 'right');
    editTime = uieditfield(pGrid, 'numeric', 'Value', 5);

    % Run Simulation Button
    btnRun = uibutton(pGrid, 'Text', 'Run', 'FontWeight', 'bold', ...
                   'ButtonPushedFcn', @(btn,event) runSim());
    btnRun.Layout.Row = 7;
    btnRun.Layout.Column = 1;

    % Stop & Close Button 
    btnStop = uibutton(pGrid, 'Text', 'Stop & Close', 'FontWeight', 'bold', ...
                   'FontColor', [0.85 0 0], ...
                   'Tooltip', 'You can also press CTRL+C', ...
                   'ButtonPushedFcn', @(btn,event) closeApp(fig, autoRunTimer));
    btnStop.Layout.Row = 7;
    btnStop.Layout.Column = 2;

    % --- Native LaTeX Equation Label ---
    % Note: The 'Interpreter' property for uilabel requires MATLAB R2022a or newer.
    eqText = '$$\dot{y} = \mathrm{SAT}(Ke) = \begin{cases} S & \mathrm{if} \; Ke \ge S \\ Ke & \mathrm{if} \; R < Ke < S \\ R & \mathrm{if} \; Ke \le R \end{cases}$$';
    eqLabel = uilabel(pGrid, 'Text', eqText, ...
                      'Interpreter', 'latex', ...
                      'HorizontalAlignment', 'center', ...
                      'VerticalAlignment', 'center');
    eqLabel.Layout.Row = 8;
    eqLabel.Layout.Column = [1 2];

    % --- Plotting Axes ---
    
    % Axes 1: Position Tracking
    ax1 = uiaxes(gl);
    ax1.Layout.Row = 1;
    ax1.Layout.Column = 2;
    title(ax1, 'Actuator Position Tracking');
    xlabel(ax1, 'Time (sec)');
    ylabel(ax1, 'Position');
    grid(ax1, 'on');

    % Axes 2: Rate Saturation
    ax2 = uiaxes(gl);
    ax2.Layout.Row = 2;
    ax2.Layout.Column = 2;
    title(ax2, 'Actuator Rate (Saturation Limits)');
    xlabel(ax2, 'Time (sec)');
    ylabel(ax2, 'Rate');
    grid(ax2, 'on');

    % 4. Start the timer now that the UI is fully built
    start(autoRunTimer);

    % --- Core Simulation Logic ---
    function runSim()
        % Fetch current parameters from UI
        K_val    = editK.Value;
        S_val    = editS.Value;
        R_val    = editR.Value;
        amp_val  = editAmp.Value;
        freq_val = editFreq.Value;
        t_end    = editTime.Value;

        % Setup ODE parameters
        tspan = [0 t_end];
        x0 = 0;
        
        % Pilot input command
        pilot_cmd = @(t) amp_val * sin(freq_val * t); 
        
        % Equations of Motion with min/max clipping
        actuator_dyn = @(t, x) max(R_val, min(S_val, K_val * (pilot_cmd(t) - x)));

        % Numerical Integration
        [t, x] = ode45(actuator_dyn, tspan, x0);

        % Post-Processing for plotting
        u = pilot_cmd(t);
        x_dot = max(R_val, min(S_val, K_val .* (u - x))); 

        % Update Position Plot (ax1)
        cla(ax1); 
        plot(ax1, t, u, '--k', 'LineWidth', 1.5);
        hold(ax1, 'on');
        plot(ax1, t, x, 'b', 'LineWidth', 1.5);
        hold(ax1, 'off');
        legend(ax1, 'Pilot Command', 'Actuator State', 'Location', 'best');

        % Update Rate Plot (ax2)
        cla(ax2);
        plot(ax2, t, x_dot, 'r', 'LineWidth', 1.5);
        hold(ax2, 'on');
        yline(ax2, S_val, '--k', 'Upper Limit');
        yline(ax2, R_val, '--k', 'Lower Limit');
        hold(ax2, 'off');
        
        % Dynamically adjust Y-limits based on active data and limits
        y_max = max(S_val, max(x_dot)) + 5;
        y_min = min(R_val, min(x_dot)) - 5;
        ylim(ax2, [y_min, y_max]);
    end

    % --- CTRL+C Keyboard Listener ---
    function handleKeyPress(fig_obj, event, timer_obj)
        % Check if the Control modifier is held down and 'c' is pressed
        if ismember('control', event.Modifier) && strcmpi(event.Key, 'c')
            closeApp(fig_obj, timer_obj);
        end
    end

    % --- Cleanup Function ---
    function closeApp(fig_obj, timer_obj)
        % Stop and delete the timer if it's valid
        if isvalid(timer_obj)
            stop(timer_obj);
            delete(timer_obj);
        end
        % Close the figure window cleanly
        delete(fig_obj);
    end
end