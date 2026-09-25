function LimitCycleAnalysisbyNIDF_GUI()
    % Create the main UI Figure
    fig = uifigure('Name', 'The Limit Cycle Analysis by Using the Negative Inverse Describing Function on the Nichols Chart', 'Position', [100, 100, 1100, 750]);
    
    % Create a main grid layout to separate the left pane (controls/info) and the right pane (plot)
    mainGrid = uigridlayout(fig, [1 2]);
    mainGrid.ColumnWidth = {360, '1x'}; 
    
    % --- Left Container (Holds Controls, Image, and Definitions) ---
    leftContainer = uigridlayout(mainGrid, [3 1]);
    leftContainer.RowHeight = {'fit', 140, '1x'}; 
    leftContainer.Padding = [0 0 0 0];
    
    % 1. Control Parameters Panel
    leftPanel = uipanel(leftContainer, 'Title', 'Control Parameters');
    controlGrid = uigridlayout(leftPanel, [10 2]); 
    controlGrid.RowHeight = repmat({'fit'}, 1, 10);
    controlGrid.ColumnWidth = {'fit', '1x'};
    controlGrid.Padding = [10 10 10 10];
    
    % Helper function to quickly create labels (with LaTeX support) and numeric fields
    function ef = createNumericField(row, latexText, htmlText, defaultVal)
        lbl = uilabel(controlGrid, 'HorizontalAlignment', 'right');
        lbl.Layout.Row = row;
        lbl.Layout.Column = 1;
        
        % If modern MATLAB version supports LaTeX in labels, use it
        if isprop(lbl, 'Interpreter')
            lbl.Interpreter = 'latex';
            lbl.Text = latexText;
            lbl.FontSize = 14;
        else
            % Fallback to HTML formatting for older MATLAB versions
            lbl.Text = htmlText;
        end
        
        ef = uieditfield(controlGrid, 'numeric', 'Value', defaultVal);
        ef.Layout.Row = row;
        ef.Layout.Column = 2;
    end
    
    % Define input fields with exact defaults and LaTeX labels
    editOmegaSp    = createNumericField(1, '$\omega_{sp}$:', '<html>&omega;<sub>sp</sub>:</html>', 2.3);
    editInvTtheta2 = createNumericField(2, '$1/T_{\theta_2}$:', '<html>1/T<sub>&theta;<sub>2</sub></sub>:</html>', 0.82);
    editZetaSp     = createNumericField(3, '$\zeta_{sp}$:', '<html>&zeta;<sub>sp</sub>:</html>', 0.3087);
    editKp         = createNumericField(4, '$K_p$:', '<html>K<sub>p</sub>:</html>', 14); 
    editMdele      = createNumericField(5, '$M_{\delta_e}$:', '<html>M<sub>&delta;<sub>e</sub></sub>:</html>', 0.537);
    editR          = createNumericField(6, 'Rate Limit ($R$):', '<html>Rate Limit (R):</html>', 15);
    editPhase      = createNumericField(7, 'Phase Shift (deg):', 'Phase Shift (deg):', -360);
    
    % Auto Zoom Checkbox
    lblZoom = uilabel(controlGrid, 'HorizontalAlignment', 'right');
    lblZoom.Layout.Row = 8;
    lblZoom.Layout.Column = 1;
    if isprop(lblZoom, 'Interpreter')
        lblZoom.Interpreter = 'latex';
        lblZoom.Text = 'Auto Zoom:';
        lblZoom.FontSize = 14;
    else
        lblZoom.Text = 'Auto Zoom:';
    end
    
    chkAutoZoom = uicheckbox(controlGrid, 'Text', '(Tight Fit to Curves)');
    chkAutoZoom.Value = 1; % On by default
    chkAutoZoom.Layout.Row = 8;
    chkAutoZoom.Layout.Column = 2;

    % Update Plot Button
    updateBtn = uibutton(controlGrid, 'Text', 'Update Plot');
    updateBtn.Layout.Row = 9;
    updateBtn.Layout.Column = [1 2];
    
    % Stop & Close Button
    closeBtn = uibutton(controlGrid, 'Text', 'Stop & Close');
    closeBtn.Layout.Row = 10;
    closeBtn.Layout.Column = [1 2];
    closeBtn.BackgroundColor = [0.85 0.3 0.3]; 
    closeBtn.FontColor = 'white';
    closeBtn.ButtonPushedFcn = @(src, event) delete(fig);
    
    % 2. Transfer Function LaTeX Panel
    formulaPanel = uipanel(leftContainer, 'Title', 'Transfer Function Formula');
    formulaGrid = uigridlayout(formulaPanel, [1 1]);
    
    tfAxes = uiaxes(formulaGrid);
    tfAxes.XColor = 'none';
    tfAxes.YColor = 'none';
    tfAxes.Color = 'none';
    tfAxes.Toolbar.Visible = 'off';
    tfAxes.Interactions = [];
    xlim(tfAxes, [0 1]);
    ylim(tfAxes, [0 1]);
    
    formulaStr = '$$\frac{\theta}{\theta_c} = Y_p(s) \frac{\theta}{\delta_e}(s) = K_p \frac{1}{s} \frac{M_{\delta_e} \left(s + \frac{1}{T_{\theta_2}}\right)}{(s^2 + 2\zeta_{sp}\omega_{sp}s + \omega_{sp}^2)}$$';
    text(tfAxes, 0.5, 0.5, formulaStr, 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 15);
    
    % 3. Parameter Definitions LaTeX Panel
    defPanel = uipanel(leftContainer, 'Title', 'Parameter Definitions');
    defGrid = uigridlayout(defPanel, [1 1]);
    
    defAxes = uiaxes(defGrid);
    defAxes.XColor = 'none';
    defAxes.YColor = 'none';
    defAxes.Color = 'none';
    defAxes.Toolbar.Visible = 'off';
    defAxes.Interactions = [];
    xlim(defAxes, [0 1]);
    ylim(defAxes, [0 1]);
    
    defText = {
        '$K_p$ : Pilot gain',
        '$R$ : Slew rate of the rate limiting element',
        '$M_{\delta_e}$ : Input signal of the rate limiting element',
        '$\omega_{sp}$ : Natural frequency',
        '$T_{\theta_2}$ : Period (sec)',
        '$\zeta_{sp}$ : Damping Ratio',
        '$\theta$ : Pitch angle',
        '$\theta_c$ : Pitch angle commanded',
        '$\delta_e$ : Elevator deflection'
    };
    
    text(defAxes, 0.02, 0.95, defText, 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'top');
    
    % --- Right Panel: UIAxes Container ---
    % Create a dedicated grid cell to hold the dynamic axes
    plotGrid = uigridlayout(mainGrid, [1 1]);
    plotGrid.Padding = [0 0 0 0];
    
    % Pack handles into a struct to pass to the update function
    appData.OmegaSp = editOmegaSp;
    appData.InvTtheta2 = editInvTtheta2;
    appData.ZetaSp = editZetaSp;
    appData.Kp = editKp;
    appData.Mdele = editMdele;
    appData.R = editR;
    appData.Phase = editPhase;
    appData.AutoZoom = chkAutoZoom;
    appData.PlotGrid = plotGrid; % Pass the container instead of the axes
    
    % Assign callback function for real-time updates and button presses
    updateFn = @(src, event) updatePlot(appData);
    
    updateBtn.ButtonPushedFcn      = updateFn;
    editOmegaSp.ValueChangedFcn    = updateFn;
    editInvTtheta2.ValueChangedFcn = updateFn;
    editZetaSp.ValueChangedFcn     = updateFn;
    editKp.ValueChangedFcn         = updateFn;
    editMdele.ValueChangedFcn      = updateFn;
    editR.ValueChangedFcn          = updateFn;
    editPhase.ValueChangedFcn      = updateFn;
    chkAutoZoom.ValueChangedFcn    = updateFn;
    
    % Call the update function once to plot the initial state
    updatePlot(appData);
end

% --- Update Logic ---
function updatePlot(appData)
    % 1. Retrieve current physical parameters from GUI
    omega_sp     = appData.OmegaSp.Value;
    inv_T_theta2 = appData.InvTtheta2.Value;
    zeta_sp      = appData.ZetaSp.Value;
    Kp           = appData.Kp.Value;
    M_del_e      = appData.Mdele.Value;
    R            = appData.R.Value;
    phaseShift   = appData.Phase.Value;
    autoZoom     = appData.AutoZoom.Value;
    
    % --- CRITICAL FIX: DESTROY AND RECREATE AXES ---
    % The Control System Toolbox attaches hidden asynchronous listeners to the axes.
    % Clearing the plot leaves these ghosts behind, which crash on the next 'drawnow'.
    % The safest workaround is to completely delete the old axes and generate a fresh one.
    delete(appData.PlotGrid.Children);
    ax = uiaxes(appData.PlotGrid);
    
    % 2. Linear system G(s) constructed dynamically from formulas
    num = Kp * M_del_e * [1, inv_T_theta2]; 
    den = [1, 2 * zeta_sp * omega_sp, omega_sp^2, 0];
    Gs = tf(num, den);
    
    % 3. Frequency range for G(jw)
    w_G = linspace(0.3, 3.5, 1000);
    
    % 4. Rate Limiter Element parameters
    Ai = Kp; 
    w_N = linspace(0.01, 10, 1001);
    w_onset = R / Ai;
    
    % 5. Describing function N(Ai,w)
    N = zeros(size(w_N));
    for k = 1:length(w_N)
        alpha = w_N(k) / w_onset;
        if alpha <= 1 % Region I
            M = 1;
            phi = 0;
        elseif alpha < 1.862 % Region II
            M = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
            phi = 0.528*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171;
        else % Region III
            varpi = w_onset / w_N(k);
            M = (4/pi)*varpi;
            phi = -acos((pi/2)*varpi);
        end
        N(k) = M * exp(1j * phi);
    end
    
    % 6. Calculate -1/N(Ai,w) and convert to FRD model
    minus_inv_N = -1 ./ N;
    response = reshape(minus_inv_N, 1, 1, []);
    sys_N = frd(response, w_N);
    
    % 7. Nichols Chart Plotting
    hold(ax, 'on');
    
    % Extract current color order so dummy legend lines match perfectly
    cOrder = colororder(ax); 
    
    % Create options object to enable the Nichols grid
    opt = nicholsoptions;
    opt.Grid = 'on';
    
    % Plot systems
    p1 = nicholsplot(ax, Gs, w_G, opt);
    p2 = nicholsplot(ax, sys_N, opt);
    
    % 8. Apply Phase Shift dynamically
    p2.PhaseMatchingEnabled = 'on';
    p2.PhaseMatchingFrequency = w_N(1);
    phase_first = rad2deg(angle(minus_inv_N(1)));
    p2.PhaseMatchingValue = phase_first + phaseShift;
    
    % Safely execute redraw now that the axes are pristine
    drawnow; 
    
    % Re-assert hold state just to ensure the axes remains stable
    hold(ax, 'on');
    
    % --- AUTO ZOOM LOGIC ---
    if autoZoom
        % Compute exact phase and magnitude bounds for G(s)
        [mag_G, phase_G] = bode(Gs, w_G);
        mag_G_dB = 20*log10(squeeze(mag_G));
        phase_G = squeeze(phase_G);
        
        % Compute exact phase and magnitude bounds for -1/N matching the chart shift
        mag_N_dB = 20*log10(abs(minus_inv_N));
        phase_N_raw = rad2deg(unwrap(angle(minus_inv_N)));
        phase_N = phase_N_raw - phase_N_raw(1) + phase_first + phaseShift;
        
        % Calculate absolute min/max to establish a tight bounding box
        min_p = min([phase_G(:); phase_N(:)]);
        max_p = max([phase_G(:); phase_N(:)]);
        min_m = min([mag_G_dB(:); mag_N_dB(:)]);
        max_m = max([mag_G_dB(:); mag_N_dB(:)]);
        
        % Add 15% visual padding to the bounding box
        dp = max(10, (max_p - min_p) * 0.15);
        dm = max(5, (max_m - min_m) * 0.15);
        
        % Apply tightly zoomed limits
        xlim(ax, [min_p - dp, max_p + dp]);
        ylim(ax, [min_m - dm, max_m + dm]);
    else
        % Let MATLAB fall back to standard wide Nichols Chart limits
        xlim(ax, 'auto');
        ylim(ax, 'auto');
    end

    % Add and capture 0 dB reference line handle
    h_yline = yline(ax, 0, '--', 'Color', [0.4 0.4 0.4]); 
    
    % Update the Title
    title(ax, 'The Limit Cycle Analysis by Using The Nichols Chart of G(j\omega) and -1/N(A_{i},\omega) of Figure 8 of Ashkenas 1964');
    
    % --- SAFE LEGEND IMPLEMENTATION ---
    h1 = plot(ax, NaN, NaN, '-', 'Color', cOrder(1, :), 'LineWidth', 1.5);
    h2 = plot(ax, NaN, NaN, '-', 'Color', cOrder(2, :), 'LineWidth', 1.5);
    legend(ax, [h1, h2, h_yline], {'G(j\omega)', '-1/N(A_i,\omega)', '0 dB'}, 'Location', 'northeast');
    
    hold(ax, 'off');
end