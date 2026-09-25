function OLOP_GUI()
    % Internal flag to handle the Stop (CTRL+C equivalent) behavior safely
    stopFlag = false;

    % Create UI Figure
    fig = uifigure('Name', 'OLOP & Jump Phenomena Analyzer (Duda''s 1997 OLOP Paper)');
    
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
    % Expanded width for left panel to 500 to ensure transfer function arrays are fully visible
    mainGrid.ColumnWidth = {500, '1x'}; 

    % Left Panel (Controls split into three sections)
    leftLayout = uigridlayout(mainGrid, [3, 1]);
    % Balanced row heights for 85% screen real estate
    leftLayout.RowHeight = {210, 180, '1x'};
    leftLayout.Padding = [0 0 0 0];

    % --- 1. Simulation Parameters Panel ---
    simPanel = uipanel(leftLayout, 'Title', 'Simulation Parameters');
    simGrid = uigridlayout(simPanel, [5, 2]);
    simGrid.ColumnWidth = {150, '1x'};
    simGrid.RowHeight = {25, 25, 25, 25, 35};

    uilabel(simGrid, 'Text', 'Pilot Cmd Amp ($q_{co}$):', 'Interpreter', 'latex', 'FontSize', 12);
    qcoEdit = uieditfield(simGrid, 'numeric', 'Value', 1.1);

    uilabel(simGrid, 'Text', 'Rate Limit ($R$):', 'Interpreter', 'latex', 'FontSize', 12);
    REdit = uieditfield(simGrid, 'numeric', 'Value', 60);

    uilabel(simGrid, 'Text', 'Gain ($K_p$):', 'Interpreter', 'latex', 'FontSize', 12);
    KpEdit = uieditfield(simGrid, 'numeric', 'Value', 13.68);

    uilabel(simGrid, 'Text', 'Resolution ($n$):', 'Interpreter', 'latex', 'FontSize', 12);
    nEdit = uieditfield(simGrid, 'numeric', 'Value', 2000, 'Limits', [100, 50000], 'RoundFractionalValues', 'on');

    % Buttons Side-by-Side
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

    % --- 2. Transfer Function Parameters Panel ---
    tfPanel = uipanel(leftLayout, 'Title', 'Transfer Function Coefficients');
    tfGrid = uigridlayout(tfPanel, [5, 2]);
    tfGrid.ColumnWidth = {80, '1x'};
    tfGrid.RowHeight = {25, 25, 25, 25, 25};

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
    
    % Split panel into 3 sections: Top JPG, Middle Eqs, Bottom PNG
    eqGrid = uigridlayout(eqPanel, [3, 1]);
    eqGrid.RowHeight = {100, '1x', 90};
    eqGrid.Padding = [5 5 5 5];

    % Top: Load the attached scheme_v0.jpg
    if isfile('scheme_v0.jpg')
        imgSchema1 = uiimage(eqGrid, 'ImageSource', 'scheme_v0.jpg', 'ScaleMethod', 'fit');
    else
        imgSchema1 = uilabel(eqGrid, 'Text', 'Save "scheme_v0.jpg" in this folder to display the block diagram.', ...
            'HorizontalAlignment', 'center', 'FontColor', [0.4 0.4 0.4], 'WordWrap', 'on');
    end
    imgSchema1.Layout.Row = 1;
    imgSchema1.Layout.Column = 1;

    % Middle: Use an invisible UI axes to render LaTeX strings reliably
    axEq = uiaxes(eqGrid);
    axEq.Layout.Row = 2;
    axEq.Layout.Column = 1;
    axEq.Visible = 'off';
    axEq.XLim = [0, 1];
    axEq.YLim = [0, 10]; 
    
    text(axEq, 0.02, 10.0, '\textbf{Controller}', 'Interpreter', 'latex', 'FontSize', 13, 'VerticalAlignment', 'top');
    txtGc  = text(axEq, 0.02, 8.6, '', 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'top');
    
    text(axEq, 0.02, 6.0, '\textbf{Rate Limiter}', 'Interpreter', 'latex', 'FontSize', 13, 'VerticalAlignment', 'top');
    strSat = '$\displaystyle \dot{y} = \mathrm{SAT}(Ke) = \left\{ \begin{array}{ll} S & \mathrm{if~} Ke \geq S \\ Ke & \mathrm{if~} R < Ke < S \\ R & \mathrm{if~} Ke \leq R \end{array} \right.$';
    text(axEq, 0.02, 4.6, strSat, 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'top');
    
    text(axEq, 0.02, 2.0, '\textbf{Longitudinal Dynamics of Aircraft}', 'Interpreter', 'latex', 'FontSize', 13, 'VerticalAlignment', 'top');
    txtGac = text(axEq, 0.02, 0.6, '', 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'top');

    % Bottom: Load the attached df.png
    if isfile('df.png')
        imgSchema2 = uiimage(eqGrid, 'ImageSource', 'df.png', 'ScaleMethod', 'fit');
    else
        imgSchema2 = uilabel(eqGrid, 'Text', 'Save "df.png" in this folder to display the describing function diagram.', ...
            'HorizontalAlignment', 'center', 'FontColor', [0.4 0.4 0.4], 'WordWrap', 'on');
    end
    imgSchema2.Layout.Row = 3;
    imgSchema2.Layout.Column = 1;

    % --- Right Panel (Plots) ---
    plotGrid = uigridlayout(mainGrid, [2, 1]);
    plotGrid.RowHeight = {'7x', '3x'}; 
    
    axNichols = uiaxes(plotGrid);
    axResidual = uiaxes(plotGrid);
    
    title(axNichols, 'Nichols Chart: Linear Loop, OLOP DF, and Stability Boundary', 'Interpreter', 'latex', 'FontSize', 13);
    title(axResidual, 'fminsearch residual vs. frequency', 'Interpreter', 'latex', 'FontSize', 13);

    % Initialize equations
    updateEquations();

    % Setup a 1-second timer to run the simulation automatically after startup
    t = timer('StartDelay', 1.0, 'TimerFcn', @(~,~) safeRunSimulation());
    start(t);

    function safeRunSimulation()
        if isvalid(fig)
            runSimulation();
        end
        if isvalid(t)
            stop(t);
            delete(t);
        end
    end

    function stopAndClose()
        stopFlag = true;
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
        statusLbl.Text = 'Status: Running (please wait)...';
        statusLbl.FontColor = [0.8 0.4 0.1];
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

            w = logspace(-1, 2, n);
            pcl = Kp*Gc / (1 + Gc*Gac);
            dtr = pi/180;

            N      = numel(w);
            deltao = zeros(1, N+1);
            phi2   = zeros(1, N+1);
            fval   = zeros(1, N);
            NoMag  = zeros(1, N);
            NoPh   = zeros(1, N);

            [magpcl, p0] = bode(pcl, w(1));
            deltao(1)    = qco * squeeze(magpcl);
            phi2(1)      = squeeze(p0);
            onset_freq   = NaN;

            for z = 1:N
                drawnow limitrate;
                if stopFlag || ~isvalid(fig)
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
            
            if stopFlag || ~isvalid(fig), return; end

            cla(axNichols);
            hold(axNichols, 'on');

            H_dol   = (10.^(NoMag/20)) .* exp(1i * NoPh);
            sys_dol = frd(H_dol(:), w, 'FrequencyUnit', 'rad/s');
            H_lin   = squeeze(freqresp(Gc*Gac, w));
            sys_lin = frd(H_lin, w, 'FrequencyUnit', 'rad/s');

            opts                      = nicholsoptions;
            opts.PhaseMatching        = 'on';
            opts.PhaseMatchingFreq    = 1;
            opts.PhaseMatchingValue   = -180;
            opts.PhaseWrapping        = 'on';
            opts.PhaseWrappingBranch  = -360;
            opts.Grid                 = 'on';

            nicholsplot(axNichols, sys_lin, 'b-', sys_dol, 'r-*', opts);

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

            v1 = [-60 -90 -100 -120 -140 -160 -180];
            v2 = [13.5 7.5 5.5 2.5 1.1 0 0];    
            plot(axNichols, v1, v2, 'k-', 'LineWidth', 2);

            h1 = plot(axNichols, NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Linear $G_c G_{ac}$');
            h2 = plot(axNichols, NaN, NaN, 'r-*', 'LineWidth', 1, 'DisplayName', 'OLOP Describing Function');
            h3 = plot(axNichols, NaN, NaN, 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm', 'DisplayName', 'OLOP DF Onset Point');
            h4 = plot(axNichols, NaN, NaN, 'k-', 'LineWidth', 2, 'DisplayName', 'Stability Boundary');
            
            legend(axNichols, [h1, h2, h3, h4], 'Location', 'best', 'Interpreter', 'latex');
            
            xlim(axNichols, [-300 -50]);
            ylim(axNichols, [-20 20]);
            grid(axNichols, 'on');
            hold(axNichols, 'off');

            cla(axResidual);
            semilogx(axResidual, w, fval, 'b-', 'LineWidth', 1.5);
            grid(axResidual, 'on');
            xlabel(axResidual, '$\omega$ [rad/s]', 'Interpreter', 'latex');
            ylabel(axResidual, 'Residual $f_{val}$', 'Interpreter', 'latex');

            if isnan(onset_freq)
                statusLbl.Text = 'Status: Complete (No onset detected)';
            else
                statusLbl.Text = sprintf('Status: Complete (Onset @ %.4g rad/s)', onset_freq);
            end
            statusLbl.FontColor = [0.1 0.6 0.1];

        catch ME
            if isvalid(fig) && ~stopFlag
                statusLbl.Text = 'Status: Error occurred!';
                statusLbl.FontColor = [1 0 0];
                uialert(fig, ME.message, 'Simulation Error');
            end
        end
    end

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