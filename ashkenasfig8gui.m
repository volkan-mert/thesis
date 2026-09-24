function ashkenasfig8gui()
    % Create the main UI Figure
    fig = uifigure('Name', 'The Negative Inverse Describing Function Analysis of Nichols Chart', 'Position', [100, 100, 1100, 750]);
    
    % Create a main grid layout to separate the left pane (controls/info) and the right pane (plot)
    mainGrid = uigridlayout(fig, [1 2]);
    mainGrid.ColumnWidth = {320, '1x'}; % 320px for the left panel container, rest for plot
    
    % --- Left Container (Holds Controls, Image, and Definitions) ---
    leftContainer = uigridlayout(mainGrid, [3 1]);
    leftContainer.RowHeight = {'fit', 120, '1x'}; % Controls fit to content, Image gets 120px, Definitions take the rest
    leftContainer.Padding = [0 0 0 0];
    
    % 1. Control Parameters Panel
    leftPanel = uipanel(leftContainer, 'Title', 'Control Parameters');
    controlGrid = uigridlayout(leftPanel, [10 2]); % Increased to 10 rows for Auto Zoom
    controlGrid.RowHeight = repmat({'fit'}, 1, 10);
    controlGrid.ColumnWidth = {'fit', '1x'};
    controlGrid.Padding = [10 10 10 10];
    
    % Helper function to quickly create labels and numeric fields
    function ef = createNumericField(row, labelText, defaultVal)
        lbl = uilabel(controlGrid, 'Text', labelText, 'HorizontalAlignment', 'right');
        lbl.Layout.Row = row;
        lbl.Layout.Column = 1;
        
        ef = uieditfield(controlGrid, 'numeric', 'Value', defaultVal);
        ef.Layout.Row = row;
        ef.Layout.Column = 2;
    end
    
    % Define input fields with default values
    editOmegaSp    = createNumericField(1, 'omega_sp:', 2.3);
    editInvTtheta2 = createNumericField(2, '1/T_theta2:', 0.82);
    editZetaSp     = createNumericField(3, 'zeta_sp:', 0.3087);
    editKp         = createNumericField(4, 'Kp:', 13.68);
    editMdele      = createNumericField(5, 'M_del_e:', 0.537);
    editR          = createNumericField(6, 'Rate Limit (R):', 15);
    editPhase      = createNumericField(7, 'Phase Shift (deg):', -360);
    
    % Auto Zoom Checkbox
    lblZoom = uilabel(controlGrid, 'Text', 'Auto Zoom:', 'HorizontalAlignment', 'right');
    lblZoom.Layout.Row = 8;
    lblZoom.Layout.Column = 1;
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
    closeBtn.BackgroundColor = [0.85 0.3 0.3]; % Red background for visibility
    closeBtn.FontColor = 'white';
    closeBtn.ButtonPushedFcn = @(src, event) delete(fig);
    
    % 2. Transfer Function Image Panel
    imgPanel = uipanel(leftContainer, 'Title', 'Transfer Function Formula');
    imgGrid = uigridlayout(imgPanel, [1 1]);
    tfImg = uiimage(imgGrid);
    tfImg.ImageSource = 'C:\Users\t0900\Documents\MATLAB\Volkan\tf.jpg'; 
    
    % 3. Parameter Definitions Panel
    defPanel = uipanel(leftContainer, 'Title', 'Parameter Definitions');
    defGrid = uigridlayout(defPanel, [1 1]);
    
    % Text array matching the specified parameter definitions exactly
    defText = {
        'Kp : pilot gain';
        'R : slew rate of the rate limiting element';
        'M_del_e : input signal of the rate limiting element';
        'omega_sp : natural frequency';
        'T_theta2 : Period (sec)';
        'zeta_sp : Damping Ratio';
        'theta : pitch angle';
        'theta_c : pitch angle commanded';
        'delta_e : elevator deflection'
    };
    
    defLabel = uilabel(defGrid, 'Text', defText);
    defLabel.VerticalAlignment = 'top';
    defLabel.WordWrap = 'on';
    
    % --- Right Panel: UIAxes ---
    ax = uiaxes(mainGrid);
    
    % Pack handles into a struct to pass to the update function
    appData.OmegaSp = editOmegaSp;
    appData.InvTtheta2 = editInvTtheta2;
    appData.ZetaSp = editZetaSp;
    appData.Kp = editKp;
    appData.Mdele = editMdele;
    appData.R = editR;
    appData.Phase = editPhase;
    appData.AutoZoom = chkAutoZoom;
    appData.Axes = ax;
    
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
    
    % Clear the axes to prevent overlapping plots
    cla(appData.Axes);
    legend(appData.Axes, 'off');
    
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
    hold(appData.Axes, 'on');
    
    % Extract current color order so dummy legend lines match perfectly
    cOrder = colororder(appData.Axes); 
    
    % Create options object to enable the Nichols grid
    opt = nicholsoptions;
    opt.Grid = 'on';
    
    p1 = nicholsplot(appData.Axes, Gs, w_G, opt);
    p2 = nicholsplot(appData.Axes, sys_N, opt);
    
    % 8. Apply Phase Shift dynamically
    p2.PhaseMatchingEnabled = 'on';
    p2.PhaseMatchingFrequency = w_N(1);
    phase_first = rad2deg(angle(minus_inv_N(1)));
    p2.PhaseMatchingValue = phase_first + phaseShift;
    
    % --- CRITICAL FIX: FORCE REDRAW ---
    drawnow; 
    
    % Re-assert hold state just to ensure the axes remains stable
    hold(appData.Axes, 'on');
    
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
        xlim(appData.Axes, [min_p - dp, max_p + dp]);
        ylim(appData.Axes, [min_m - dm, max_m + dm]);
    else
        % Let MATLAB fall back to standard wide Nichols Chart limits
        xlim(appData.Axes, 'auto');
        ylim(appData.Axes, 'auto');
    end

    % Add and capture 0 dB reference line handle
    h_yline = yline(appData.Axes, 0, '--', 'Color', [0.4 0.4 0.4]); 
    
    % Update the Title
    title(appData.Axes, 'The Nichols Chart of G(j\omega) and -1/N(A_{i},\omega) of Figure 8 of Ashkenas 1964');
    
    % --- SAFE LEGEND IMPLEMENTATION ---
    h1 = plot(appData.Axes, NaN, NaN, '-', 'Color', cOrder(1, :), 'LineWidth', 1.5);
    h2 = plot(appData.Axes, NaN, NaN, '-', 'Color', cOrder(2, :), 'LineWidth', 1.5);
    legend(appData.Axes, [h1, h2, h_yline], {'G(j\omega)', '-1/N(A_i,\omega)', '0 dB'}, 'Location', 'northeast');
    
    hold(appData.Axes, 'off');
end