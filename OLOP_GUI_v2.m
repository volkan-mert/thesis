function OLOP_GUI_v2()
    % Internal flags to handle safe execution and stopping
    stopFlag = false;
    isComputing = false; % Prevents concurrent overlapping executions

    % Create UI Figure
    fig = uifigure('Name', 'Comprehensive OLOP & PIO Analyzer');
    
    % Dynamically calculate 85% of the screen size and center it
    screenSize = get(groot, 'ScreenSize');
    figWidth  = screenSize(3) * 0.85;
    figHeight = screenSize(4) * 0.85;
    figX      = (screenSize(3) - figWidth) / 2;
    figY      = (screenSize(4) - figHeight) / 2;
    
    % Apply the calculated position and size
    fig.Position = [figX, figY, figWidth, figHeight];
    
    % Ensure the figure handles standard closing gracefully
    fig.CloseRequestFcn = @(src, event) stopAndClose();

    % Main Layout
    mainGrid = uigridlayout(fig, [1, 2]);
    mainGrid.ColumnWidth = {500, '1x'}; 

    % Left Panel
    leftLayout = uigridlayout(mainGrid, [3, 1]);
    leftLayout.RowHeight = {260, 250, '1x'};
    leftLayout.Padding = [0 0 0 0];

    % --- 1. Simulation Parameters Panel ---
    simPanel = uipanel(leftLayout, 'Title', 'Simulation Parameters');
    simGrid = uigridlayout(simPanel, [6, 2]); 
    simGrid.ColumnWidth = {150, '1x'};
    simGrid.RowHeight = {25, 25, 25, 25, 40, 30};

    uilabel(simGrid, 'Text', 'Pilot Cmd Amp ($q_{co}$):', 'Interpreter', 'latex', 'FontSize', 12);
    qcoEdit = uieditfield(simGrid, 'numeric', 'Value', 1.1);

    uilabel(simGrid, 'Text', 'Rate Limit ($R$):', 'Interpreter', 'latex', 'FontSize', 12);
    REdit = uieditfield(simGrid, 'numeric', 'Value', 60);

    uilabel(simGrid, 'Text', 'Gain ($K_p$):', 'Interpreter', 'latex', 'FontSize', 12);
    KpEdit = uieditfield(simGrid, 'numeric', 'Value', 13.68);

    uilabel(simGrid, 'Text', 'Resolution ($n$):', 'Interpreter', 'latex', 'FontSize', 12);
    nEdit = uieditfield(simGrid, 'numeric', 'Value', 2000, 'Limits', [100, 50000], 'RoundFractionalValues', 'on');

    runBtn = uibutton(simGrid, 'Text', 'Run Analysis', 'ButtonPushedFcn', @(src, event) runSimulation());
    runBtn.Layout.Column = 1;
    runBtn.BackgroundColor = [0 0.447 0.741]; 
    runBtn.FontColor = [1 1 1];               
    runBtn.FontWeight = 'bold';
    
    stopBtn = uibutton(simGrid, 'Text', 'Stop & Close', 'ButtonPushedFcn', @(src, event) stopAndClose());
    stopBtn.Layout.Column = 2;
    stopBtn.BackgroundColor = [0.8 0.2 0.2]; 
    stopBtn.FontColor = [1 1 1];             
    stopBtn.FontWeight = 'bold';

    % Automated Verdict Indicator
    vGrid = uigridlayout(simGrid, [1 2]);
    vGrid.Layout.Row = 6;
    vGrid.Layout.Column = [1 2];
    vGrid.ColumnWidth = {30, '1x'};
    vGrid.Padding = [0 0 0 0];
    
    verdictLamp = uilamp(vGrid);
    verdictLamp.Color = [0.5 0.5 0.5];
    verdictLbl = uilabel(vGrid, 'Text', 'PIO Verdict: Standby...', 'FontWeight', 'bold', 'FontSize', 13);

    % --- 2. Transfer Function Parameters Panel ---
    tfPanel = uipanel(leftLayout, 'Title', 'Transfer Function Coefficients');
    tfGrid = uigridlayout(tfPanel, [5, 2]);
    tfGrid.ColumnWidth = {110, '1x'};
    tfGrid.RowHeight = {30, 30, 30, 30, 30};

    uilabel(tfGrid, 'Text', '$\mathrm{num}(G_c)$:', 'Interpreter', 'latex', 'FontSize', 13);
    numGcEdit = uieditfield(tfGrid, 'text', 'Value', '[5.21, -273.7855, -1425.240306, -700.1952408]');
    
    uilabel(tfGrid, 'Text', '$\mathrm{den}(G_c)$:', 'Interpreter', 'latex', 'FontSize', 13);
    denGcEdit = uieditfield(tfGrid, 'text', 'Value', '[1, 21.3594, 545.553804, 605.6621, 0]');

    uilabel(tfGrid, 'Text', '$\mathrm{num}(G_{ac})$:', 'Interpreter', 'latex', 'FontSize', 13);
    numGacEdit = uieditfield(tfGrid, 'text', 'Value', '[-10.524, -16.8384, -0.62466254, 0]');

    uilabel(tfGrid, 'Text', '$\mathrm{den}(G_{ac})$:', 'Interpreter', 'latex', 'FontSize', 13);
    denGacEdit = uieditfield(tfGrid, 'text', 'Value', '[1, 2.347312, -5.30606528, -0.18359616, -0.0418176]');

    statusLbl = uilabel(tfGrid, 'Text', 'Status: Ready');
    statusLbl.Layout.Column = [1 2];
    statusLbl.FontWeight = 'bold';

    numGcEdit.ValueChangedFcn = @(src, event) updateEquations();
    denGcEdit.ValueChangedFcn = @(src, event) updateEquations();
    numGacEdit.ValueChangedFcn = @(src, event) updateEquations();
    denGacEdit.ValueChangedFcn = @(src, event) updateEquations();

    % --- 3. Equations & Block Diagram Display Panel ---
    eqPanel = uipanel(leftLayout, 'Title', 'System Equations');
    eqGrid = uigridlayout(eqPanel, [3, 1]);
    eqGrid.RowHeight = {100, '1x', 90};
    eqGrid.Padding = [5 5 5 5];

    if isfile('scheme_v0.jpg')
        imgSchema1 = uiimage(eqGrid, 'ImageSource', 'scheme_v0.jpg', 'ScaleMethod', 'fit');
    else
        imgSchema1 = uilabel(eqGrid, 'Text', 'Save "scheme_v0.jpg" in this folder.', 'HorizontalAlignment', 'center');
    end
    imgSchema1.Layout.Row = 1;
    imgSchema1.Layout.Column = 1;

    axEq = uiaxes(eqGrid);
    axEq.Layout.Row = 2;
    axEq.Layout.Column = 1;
    axEq.Visible = 'off';
    axEq.XLim = [0, 1];
    axEq.YLim = [0, 10]; 
    
    % PREVENTS THE UI INTERACTION CRASH ON INVISIBLE AXES
    axEq.HitTest = 'off';
    axEq.PickableParts = 'none';
    
    text(axEq, 0.02, 10.0, '\textbf{Controller}', 'Interpreter', 'latex', 'FontSize', 13, 'VerticalAlignment', 'top');
    txtGc  = text(axEq, 0.02, 8.6, '', 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'top');
    
    text(axEq, 0.02, 6.0, '\textbf{Rate Limiter}', 'Interpreter', 'latex', 'FontSize', 13, 'VerticalAlignment', 'top');
    strSat = '$\displaystyle \dot{y} = \mathrm{SAT}(Ke) = \left\{ \begin{array}{ll} S & \mathrm{if~} Ke \geq S \\ Ke & \mathrm{if~} R < Ke < S \\ R & \mathrm{if~} Ke \leq R \end{array} \right.$';
    text(axEq, 0.02, 4.6, strSat, 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'top');
    
    text(axEq, 0.02, 2.0, '\textbf{Longitudinal Dynamics of Aircraft}', 'Interpreter', 'latex', 'FontSize', 13, 'VerticalAlignment', 'top');
    txtGac = text(axEq, 0.02, 0.6, '', 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'top');

    if isfile('df.png')
        imgSchema2 = uiimage(eqGrid, 'ImageSource', 'df.png', 'ScaleMethod', 'fit');
    else
        imgSchema2 = uilabel(eqGrid, 'Text', 'Save "df.png" in this folder.', 'HorizontalAlignment', 'center');
    end
    imgSchema2.Layout.Row = 3;
    imgSchema2.Layout.Column = 1;

    % --- Right Panel (Tabs for Multiple Diagnostics) ---
    tabGroup = uitabgroup(mainGrid);
    
    % TAB 1: Nichols Plot
    tabNichols = uitab(tabGroup, 'Title', 'OLOP Nichols Chart');
    gridNichols = uigridlayout(tabNichols, [2 1]);
    gridNichols.RowHeight = {'8x', '2x'}; 
    axNichols = uiaxes(gridNichols);
    axResidual = uiaxes(gridNichols);
    
    % TAB 2: Actuator State
    tabActuator = uitab(tabGroup, 'Title', 'Actuator Saturation');
    gridAct = uigridlayout(tabActuator, [2 1]);
    axDelta = uiaxes(gridAct);
    axXparam = uiaxes(gridAct);
    
    % TAB 3: Describing Function
    tabDF = uitab(tabGroup, 'Title', 'Describing Function');
    gridDF = uigridlayout(tabDF, [2 1]);
    axDFMag = uiaxes(gridDF);
    axDFPh = uiaxes(gridDF);
    
    % TAB 4: Linear Bode
    tabBode = uitab(tabGroup, 'Title', 'Linear Baseline (Bode)');
    gridBode = uigridlayout(tabBode, [2 1]);
    axBodeMag = uiaxes(gridBode);
    axBodePh = uiaxes(gridBode);

    % Formats
    title(axNichols, 'Nichols Chart', 'Interpreter', 'latex', 'FontSize', 13);
    title(axResidual, 'Harmonic Balance Residual vs. Frequency', 'Interpreter', 'latex', 'FontSize', 13);
    
    title(axDelta, 'Actuator Surface-Command Amplitude ($\delta_0$)', 'Interpreter', 'latex', 'FontSize', 13);
    title(axXparam, 'Saturation State Parameter ($x = \omega \delta_0 / R$)', 'Interpreter', 'latex', 'FontSize', 13);
    
    title(axDFMag, 'Describing Function Gain ($N_{mag}$)', 'Interpreter', 'latex', 'FontSize', 13);
    title(axDFPh, 'Describing Function Phase ($N_{ph}$)', 'Interpreter', 'latex', 'FontSize', 13);
    
    title(axBodeMag, 'Linear System Baseline: Magnitude', 'Interpreter', 'latex', 'FontSize', 13);
    title(axBodePh, 'Linear System Baseline: Phase', 'Interpreter', 'latex', 'FontSize', 13);

    % Initialize equations and auto-start
    updateEquations();
    t = timer('StartDelay', 1.0, 'TimerFcn', @(~,~) safeRunSimulation());
    start(t);

    function safeRunSimulation()
        if isvalid(fig) && ~isComputing
            runSimulation();
        end
        if isvalid(t)
            stop(t);
            delete(t);
        end
    end

    function stopAndClose()
        stopFlag = true;
        if isvalid(t)
            stop(t);
            delete(t);
        end
        if isvalid(fig)
            delete(fig);
        end
    end

    function updateEquations()
        num_Gc_val  = str2num(numGcEdit.Value); %#ok<ST2NM>
        den_Gc_val  = str2num(denGcEdit.Value); %#ok<ST2NM>
        num_Gac_val = str2num(numGacEdit.Value); %#ok<ST2NM>
        den_Gac_val = str2num(denGacEdit.Value); %#ok<ST2NM>
        
        if isempty(num_Gc_val), num_Gc_val = 0; end
        if isempty(den_Gc_val), den_Gc_val = 1; end
        if isempty(num_Gac_val), num_Gac_val = 0; end
        if isempty(den_Gac_val), den_Gac_val = 1; end

        strGc  = sprintf('$\\displaystyle G_c(s) = \\frac{%s}{%s} $', poly2latex(num_Gc_val), poly2latex(den_Gc_val));
        strGac = sprintf('$\\displaystyle G_{ac}(s) = \\frac{%s}{%s} $', poly2latex(num_Gac_val), poly2latex(den_Gac_val));
        
        txtGc.String  = strGc;
        txtGac.String = strGac;
    end

    function runSimulation()
        % Prevent overlapping executions
        if isComputing || ~isvalid(fig)
            return;
        end
        isComputing = true;
        
        statusLbl.Text = 'Status: Running (please wait)...';
        statusLbl.FontColor = [0.8 0.4 0.1];
        verdictLamp.Color = [0.5 0.5 0.5];
        verdictLbl.Text = 'PIO Verdict: Calculating...';
        verdictLbl.FontColor = [0 0 0];
        stopFlag = false; 
        drawnow; 

        try
            qco = qcoEdit.Value;
            R   = REdit.Value;
            Kp  = KpEdit.Value;
            n   = nEdit.Value;

            num_Gc_val  = str2num(numGcEdit.Value); %#ok<ST2NM>
            den_Gc_val  = str2num(denGcEdit.Value); %#ok<ST2NM>
            num_Gac_val = str2num(numGacEdit.Value); %#ok<ST2NM>
            den_Gac_val = str2num(denGacEdit.Value); %#ok<ST2NM>

            Gc = tf(num_Gc_val, den_Gc_val);
            Gac = tf(num_Gac_val, den_Gac_val);

            % Linear Margins Calculation
            [Gm, Pm, Wcg, Wcp] = margin(Gc*Gac);
            Gm_dB = 20*log10(Gm);
            if isinf(Gm_dB)
                gm_str = '\infty';
            else
                gm_str = sprintf('%.2f', Gm_dB);
            end

            w = logspace(-1, 2, n);
            pcl = Kp*Gc / (1 + Gc*Gac);
            dtr = pi/180;

            N        = numel(w);
            deltao   = zeros(1, N+1);
            phi2     = zeros(1, N+1);
            fval     = zeros(1, N);
            NoMag    = zeros(1, N);
            NoPh     = zeros(1, N);
            
            % Arrays for new plots
            magN_arr = zeros(1, N);
            phN_arr  = zeros(1, N);
            x_arr    = zeros(1, N);

            [magpcl, p0] = bode(pcl, w(1));
            deltao(1)    = qco * squeeze(magpcl);
            phi2(1)      = squeeze(p0);
            onset_freq   = NaN;

            for z = 1:N
                drawnow limitrate;
                if stopFlag || ~isvalid(fig)
                    isComputing = false;
                    return; 
                end

                [mGc,  pGc ] = bode(Gc,  w(z));
                [mGac, pGac] = bode(Gac, w(z));
                magGc  = squeeze(mGc);   phGc  = squeeze(pGc);
                magGac = squeeze(mGac);  phGac = squeeze(pGac);

                xo = [deltao(z), phi2(z)];
                
                obj_fun = @(x) eqs(x, w(z), magGc, phGc, magGac, phGac, Kp, qco, R);
                [a, fval(z)] = fminsearch(obj_fun, xo, optimset('Display','off'));

                deltait = a(1);
                phi2it  = a(2);
                [magNit, phNit] = dfunction(w(z), R, deltait);
                
                % Store descriptive data
                magN_arr(z) = magNit;
                phN_arr(z)  = phNit;
                x_arr(z)    = w(z) * deltait / R;
                
                A   = deltait * magGac * magNit / (Kp * qco);
                phi = phi2it + phGac + phNit*180/pi;
                
                NoMag(z) = 20*log10( A / sqrt(1 - 2*A*cos(dtr*phi) + A^2) );
                NoPh(z)  = dtr*phi - atan2(-A*sin(dtr*phi), 1 - A*cos(dtr*phi));

                if fval(z) > 1e-4
                    if isnan(onset_freq)
                        onset_freq = w(z);
                    end
                    deltao(z+1) = 24;       
                    phi2(z+1)   = -237;   
                else
                    deltao(z+1) = a(1);
                    phi2(z+1)   = a(2);
                end
            end
            
            if stopFlag || ~isvalid(fig)
                isComputing = false;
                return; 
            end

            % -------------------------------------------------------------
            % TAB 1: Nichols Plot Updating
            % -------------------------------------------------------------
            cla(axNichols);
            H_dol   = (10.^(NoMag/20)) .* exp(1i * NoPh);
            sys_dol = frd(H_dol(:), w, 'FrequencyUnit', 'rad/s');
            H_lin   = squeeze(freqresp(Gc*Gac, w));
            sys_lin = frd(H_lin, w, 'FrequencyUnit', 'rad/s');

            opts = nicholsoptions;
            opts.PhaseMatching       = 'on';
            opts.PhaseMatchingFreq   = 1;
            opts.PhaseMatchingValue  = -180;
            opts.PhaseWrapping       = 'on';
            opts.PhaseWrappingBranch = -360;
            opts.Grid                = 'on';

            nicholsplot(axNichols, sys_lin, 'b-', sys_dol, 'r-*', opts);
            hold(axNichols, 'on');
            
            % Stability Boundary Vectors
            v1 = [-60 -90 -100 -120 -140 -160 -180];
            v2 = [13.5 7.5 5.5 2.5 1.1 0 0]; 

            mag_onset_dB = NaN;
            ph_onset_deg = NaN;

            if ~isnan(onset_freq)
                resp = squeeze(freqresp(sys_dol, onset_freq));
                mag_onset_dB = 20*log10(abs(resp));
                ph_onset_deg = rad2deg(angle(resp));
                
                while ph_onset_deg > 0
                    ph_onset_deg = ph_onset_deg - 360;
                end
                while ph_onset_deg < -360
                    ph_onset_deg = ph_onset_deg + 360;
                end
                
                plot(axNichols, ph_onset_deg, mag_onset_dB, 'mp', 'MarkerSize', 16, 'MarkerFaceColor', 'm', 'HandleVisibility', 'off');
            end

            xline(axNichols, -180, 'k--', 'HandleVisibility', 'off');
            yline(axNichols, 0, '--', 'HandleVisibility', 'off');
            plot(axNichols, -180, 0, 'r+', 'HandleVisibility', 'off');   
            
            % Plot the Stability Boundary
            plot(axNichols, v1, v2, 'k-', 'LineWidth', 2);

            % --- Annotate Stable and Unstable Regions ---
            text(axNichols, -110, 8, '\textbf{Unstable Region}', 'Interpreter', 'latex', 'FontSize', 14, 'Color', [0.7 0 0], 'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
            text(axNichols, -130, -5, '\textbf{Stable Region}', 'Interpreter', 'latex', 'FontSize', 14, 'Color', [0 0.5 0], 'HorizontalAlignment', 'center', 'HandleVisibility', 'off');

            % Update Titles and overlay Margins
            title(axNichols, sprintf('Nichols Chart (GM: %s dB, PM: %.1f$^{\\circ}$)', gm_str, Pm), 'Interpreter', 'latex', 'FontSize', 13);
            
            h1 = plot(axNichols, NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Linear $G_c G_{ac}$');
            h2 = plot(axNichols, NaN, NaN, 'r-*', 'LineWidth', 1, 'DisplayName', 'OLOP Describing Function');
            h3 = plot(axNichols, NaN, NaN, 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm', 'DisplayName', 'OLOP DF Onset Point');
            h4 = plot(axNichols, NaN, NaN, 'k-', 'LineWidth', 2, 'DisplayName', 'Stability Boundary');
            
            % Plot GM and PM physical marks on Nichols
            h5 = plot(axNichols, NaN, NaN, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', sprintf('Gain Margin (%s dB)', gm_str));
            h6 = plot(axNichols, NaN, NaN, 'gd', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', sprintf('Phase Margin (%.1f$^{\\circ}$)', Pm));
            
            if ~isinf(Gm_dB) && ~isnan(Gm_dB) && ~isnan(Wcg)
                plot(axNichols, -180, -Gm_dB, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
            end
            if ~isnan(Pm) && ~isnan(Wcp)
                plot(axNichols, -180 + Pm, 0, 'gd', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
            end

            legend(axNichols, [h1, h2, h3, h4, h5, h6], 'Location', 'best', 'Interpreter', 'latex');
            xlim(axNichols, [-300 -50]);
            ylim(axNichols, [-20 20]);
            grid(axNichols, 'on');
            hold(axNichols, 'off');

            cla(axResidual);
            semilogx(axResidual, w, fval, 'b-', 'LineWidth', 1.5);
            grid(axResidual, 'on');
            xlabel(axResidual, '$\omega$ [rad/s]', 'Interpreter', 'latex');
            ylabel(axResidual, 'Residual $f_{val}$', 'Interpreter', 'latex');

            % --- PIO Verdict Logic ---
            if isnan(onset_freq)
                verdictLamp.Color = [0 0.8 0]; % Green
                verdictLbl.Text = 'SAFE: No Saturation / No Onset Point Found';
                verdictLbl.FontColor = [0 0.6 0];
            else
                if ph_onset_deg >= -180 && ph_onset_deg <= -60
                    bound_mag = interp1(v1, v2, ph_onset_deg, 'linear', 'extrap');
                    if mag_onset_dB > bound_mag
                        verdictLamp.Color = [0.8 0 0]; % Red
                        verdictLbl.Text = 'WARNING: Category II PIO Predicted!';
                        verdictLbl.FontColor = [0.8 0 0];
                    else
                        verdictLamp.Color = [0 0.8 0]; % Green
                        verdictLbl.Text = 'SAFE: Onset lies below Stability Boundary';
                        verdictLbl.FontColor = [0 0.6 0];
                    end
                elseif ph_onset_deg < -180
                    verdictLamp.Color = [0.8 0 0]; % Red
                    verdictLbl.Text = 'WARNING: Phase < -180 at Onset. PIO Likely!';
                    verdictLbl.FontColor = [0.8 0 0];
                else
                    verdictLamp.Color = [0 0.8 0]; % Green
                    verdictLbl.Text = 'SAFE: Onset point has adequate phase margin';
                    verdictLbl.FontColor = [0 0.6 0];
                end
            end

            % -------------------------------------------------------------
            % TAB 2: Actuator State updating
            % -------------------------------------------------------------
            cla(axDelta);
            loglog(axDelta, w, deltao(1:N), 'b-', 'LineWidth', 1.5);
            hold(axDelta, 'on');
            loglog(axDelta, w, R./w, 'r--', 'LineWidth', 1.5);
            grid(axDelta, 'on');
            legend(axDelta, 'Cmd Amp ($\delta_0$)', 'Saturation Bound ($R/\omega$)', 'Location', 'best', 'Interpreter', 'latex');
            xlabel(axDelta, '$\omega$ [rad/s]', 'Interpreter', 'latex');
            ylabel(axDelta, 'Amplitude [deg]', 'Interpreter', 'latex');
            hold(axDelta, 'off');

            cla(axXparam);
            semilogx(axXparam, w, x_arr, 'k-', 'LineWidth', 1.5);
            hold(axXparam, 'on');
            yline(axXparam, 1, 'g--', 'Linear Boundary ($x=1$)', 'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left');
            yline(axXparam, 1.862, 'r--', 'Fully Saturated ($x=1.862$)', 'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left');
            grid(axXparam, 'on');
            xlabel(axXparam, '$\omega$ [rad/s]', 'Interpreter', 'latex');
            ylabel(axXparam, '$x$', 'Interpreter', 'latex');
            hold(axXparam, 'off');

            % -------------------------------------------------------------
            % TAB 3: Describing Function updating
            % -------------------------------------------------------------
            cla(axDFMag);
            semilogx(axDFMag, w, magN_arr, 'b-', 'LineWidth', 1.5);
            grid(axDFMag, 'on');
            xlabel(axDFMag, '$\omega$ [rad/s]', 'Interpreter', 'latex');
            ylabel(axDFMag, 'Gain', 'Interpreter', 'latex');

            cla(axDFPh);
            semilogx(axDFPh, w, rad2deg(phN_arr), 'r-', 'LineWidth', 1.5);
            grid(axDFPh, 'on');
            xlabel(axDFPh, '$\omega$ [rad/s]', 'Interpreter', 'latex');
            ylabel(axDFPh, 'Phase [deg]', 'Interpreter', 'latex');

            % -------------------------------------------------------------
            % TAB 4: Linear Bode updating
            % -------------------------------------------------------------
            [magOL, phOL] = bode(Gc*Gac, w);
            [magCL, phCL] = bode(pcl, w);
            
            title(axBodeMag, sprintf('Linear System Baseline: Magnitude (GM: %s dB, PM: %.1f$^{\\circ}$)', gm_str, Pm), 'Interpreter', 'latex', 'FontSize', 13);
            cla(axBodeMag);
            semilogx(axBodeMag, w, 20*log10(squeeze(magOL)), 'b-', 'LineWidth', 1.5);
            hold(axBodeMag, 'on');
            semilogx(axBodeMag, w, 20*log10(squeeze(magCL)), 'm--', 'LineWidth', 1.5);
            
            % Overlay GM and PM markers and lines on Bode Magnitude
            if ~isinf(Gm_dB) && ~isnan(Gm_dB) && ~isnan(Wcg)
                xline(axBodeMag, Wcg, 'r:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                plot(axBodeMag, Wcg, -Gm_dB, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
            end
            if ~isnan(Pm) && ~isnan(Wcp)
                xline(axBodeMag, Wcp, 'g:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                plot(axBodeMag, Wcp, 0, 'gd', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
            end

            grid(axBodeMag, 'on');
            legend(axBodeMag, 'Open Loop ($G_c G_{ac}$)', 'Closed Loop', 'Interpreter', 'latex', 'Location', 'best');
            xlabel(axBodeMag, '$\omega$ [rad/s]', 'Interpreter', 'latex');
            ylabel(axBodeMag, 'Magnitude [dB]', 'Interpreter', 'latex');
            hold(axBodeMag, 'off');

            cla(axBodePh);
            semilogx(axBodePh, w, squeeze(phOL), 'b-', 'LineWidth', 1.5);
            hold(axBodePh, 'on');
            semilogx(axBodePh, w, squeeze(phCL), 'm--', 'LineWidth', 1.5);
            
            % Overlay GM and PM markers and lines on Bode Phase
            if ~isinf(Gm_dB) && ~isnan(Gm_dB) && ~isnan(Wcg)
                xline(axBodePh, Wcg, 'r:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                plot(axBodePh, Wcg, -180, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
            end
            if ~isnan(Pm) && ~isnan(Wcp)
                xline(axBodePh, Wcp, 'g:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                plot(axBodePh, Wcp, -180 + Pm, 'gd', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
            end

            grid(axBodePh, 'on');
            xlabel(axBodePh, '$\omega$ [rad/s]', 'Interpreter', 'latex');
            ylabel(axBodePh, 'Phase [deg]', 'Interpreter', 'latex');
            hold(axBodePh, 'off');

            statusLbl.Text = 'Status: Analysis Complete';
            statusLbl.FontColor = [0.1 0.6 0.1];
            
            isComputing = false; % Release lock

        catch ME
            isComputing = false; % Release lock on error
            if isvalid(fig) && ~stopFlag
                statusLbl.Text = 'Status: Error occurred!';
                statusLbl.FontColor = [1 0 0];
                uialert(fig, ME.message, 'Simulation Error');
            end
        end
    end

    % --- Helper Math Functions ---
    function str = poly2latex(coeffs)
        if isempty(coeffs)
            str = '0'; return;
        end
        str = '';
        degree = length(coeffs) - 1;
        
        for i = 1:length(coeffs)
            c = coeffs(i);
            if c == 0 && degree > 0
                continue;
            end
            
            if c > 0 && ~isempty(str)
                str = [str ' + ']; %#ok<AGROW>
            elseif c < 0
                if isempty(str)
                    str = '-';
                else
                    str = [str ' - ']; %#ok<AGROW>
                end
            end
            
            abs_c = abs(c);
            if abs_c ~= 1 || (degree - i + 1) == 0
                str = [str num2str(abs_c, '%.4g')]; %#ok<AGROW>
            end
            
            p = degree - i + 1;
            if p > 1
                str = [str 's^{' num2str(p) '}']; %#ok<AGROW>
            elseif p == 1
                str = [str 's']; %#ok<AGROW>
            end
        end
        if isempty(str), str = '0'; end
    end

    function f = eqs(x, freq, mag_Gc, ph_Gc, mag_Gac, ph_Gac, Kp_val, qco_val, R_val)
        deltao = x(1);
        phi2   = x(2);
        dtr    = pi/180;
        
        [magN, phN] = dfunction(freq, R_val, deltao);
        
        re = deltao/mag_Gc * cos(dtr*(phi2 - ph_Gc)) + deltao*mag_Gac*magN * cos(dtr*(phi2 + ph_Gac) + phN) - Kp_val*qco_val;
        im = deltao/mag_Gc * sin(dtr*(phi2 - ph_Gc)) + deltao*mag_Gac*magN * sin(dtr*(phi2 + ph_Gac) + phN);
        f  = re^2 + im^2;
    end

    function [magN, phN] = dfunction(freq, rate, inamp)
        x_val = freq * inamp / rate;
        if x_val < 1
            magN = 1;
            phN  = 0;
        elseif x_val < 1.862
            magN = polyval([0.2908 -1.4396 1.9232 0.2230], x_val);
            phN  = polyval([0.5280 -2.6213 3.5056 -1.4171], x_val);
        else
            magN = 4 / (x_val*pi);
            phN  = -acos(pi/(2*x_val));
        end
    end
end