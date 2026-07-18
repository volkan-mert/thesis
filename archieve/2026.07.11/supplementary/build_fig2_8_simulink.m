%% ========================================================================
%  build_fig2_8_simulink.m
%
%  Programmatically builds a Simulink model that reproduces Figure 2-8
%  "Time Domain Analysis" from Gilbreath (2001), OLOP / actuator-rate-limit
%  PIO thesis, then runs it at omega = 5.0 and 5.3 rad/sec and plots the
%  result in the same 2x2 layout as the thesis figure.
%
%  Closed loop (Fig. 2-3, pilot loop open):
%
%     qc --[K]--(+)-->[ Gc ]-->[ Rate Limiter (+/-R) ]-->[ Gac ]--+--> q
%                (-)                                               |
%                 ^________________________________________________|
%
%  Model data are verbatim from the thesis Appendix B M-file
%  (AIAA-95-3204-CP example).  Onset frequency ~5.1 rad/sec:
%     omega = 5.0  -> bounded / stable
%     omega = 5.3  -> rate limiter activates, loop diverges after ~3 s
%
%  Requires: Simulink + Control System Toolbox.  Run as a SCRIPT.
%% ========================================================================
clear; clc; close all; bdclose all;

%% ---- 1. System data (base workspace vars referenced by the TF blocks) ---
K   = 13.68;      % command / pilot-station gain
R   = 60;         % actuator rate limit                        [deg/sec]
qco = 1.1;        % pitch-rate command amplitude               [deg/sec]

numGc  = 5.21 * conv(conv([1 -57.36],[1 4.26]),[1 0.55]);
denGc  = conv(conv([1 2*0.442*22.85 22.85^2],[1 0]),[1 1.16]);

numGac = -10.524 * conv(conv([1 1.562],[1 0.038]),[1 0]);
denGac = conv(conv([1 2*0.212*0.088 0.088^2],[1 3.75]),[1 -1.44]);

%% ---- 2. Build the model -------------------------------------------------
mdl = 'fig2_8_pio';
if bdIsLoaded(mdl); close_system(mdl,0); end
new_system(mdl);  open_system(mdl);

add = @(src,name,pos,varargin) add_block(src,[mdl '/' name], ...
        'Position',pos,varargin{:});

% Blocks  (positions: [left top right bottom])
add('simulink/Sources/Sine Wave','qc',[30 150 70 190], ...
        'Amplitude','qco','Frequency','w','Phase','0');
add('simulink/Math Operations/Gain','K',[110 152 150 188],'Gain','K');
add('simulink/Math Operations/Sum','Sum',[190 152 220 188], ...
        'Inputs','+-','IconShape','round');
add('simulink/Continuous/Transfer Fcn','Gc',[260 145 360 195], ...
        'Numerator','numGc','Denominator','denGc');
add('simulink/Discontinuities/Rate Limiter','RLE',[400 148 450 192], ...
        'RisingSlewLimit','R','FallingSlewLimit','-R');
add('simulink/Continuous/Transfer Fcn','Gac',[490 145 590 195], ...
        'Numerator','numGac','Denominator','denGac');

% Logging (To Workspace, timeseries)
add('simulink/Sinks/To Workspace','q_log' ,[650 150 710 190], ...
        'VariableName','q_log' ,'SaveFormat','Timeseries');
add('simulink/Sinks/To Workspace','cmd_log',[190 60 250 100], ...
        'VariableName','cmd_log','SaveFormat','Timeseries');   % = K*qc
add('simulink/Sinks/To Workspace','dc_log',[400 60 460 100], ...
        'VariableName','dc_log','SaveFormat','Timeseries');    % delta_c
add('simulink/Sinks/To Workspace','dlt_log',[490 250 550 290], ...
        'VariableName','dlt_log','SaveFormat','Timeseries');   % delta

% Wiring (autorouting handles the branches/feedback)
L = @(a,b) add_line(mdl,a,b,'autorouting','on');
L('qc/1','K/1');
L('K/1','Sum/1');          L('K/1','cmd_log/1');   % branch: command ref
L('Sum/1','Gc/1');
L('Gc/1','RLE/1');         L('Gc/1','dc_log/1');   % branch: delta_c
L('RLE/1','Gac/1');        L('RLE/1','dlt_log/1'); % branch: delta
L('Gac/1','q_log/1');
L('Gac/1','Sum/2');                                % unity feedback

%% ---- 3. Solver configuration -------------------------------------------
% Stiff, variable-step: Gc has poles near 22.85 rad/s and an RHP zero,
% Gac has an unstable pole at +1.44. Small MaxStep resolves the limiter.
set_param(mdl, ...
    'Solver','ode15s', 'StopTime','5', ...
    'MaxStep','0.005', 'RelTol','1e-6', 'AbsTol','1e-8', ...
    'ReturnWorkspaceOutputs','on');

%% ---- 4. Run both cases and plot ----------------------------------------
omegas = [5.0 5.3];

fig = figure('Color','w','Position',[80 80 1050 680]);
tl  = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:numel(omegas)
    w = omegas(j);                 %#ok<NASGU>  (referenced by Sine block)
    assignin('base','w',w);

    simOut = sim(mdl);

    q   = simOut.q_log;   cmd = simOut.cmd_log;
    dc  = simOut.dc_log;  dlt = simOut.dlt_log;

    % ---- top row: pitch rate ----
    nexttile(j);
    plot(q.Time, q.Data, 'b-','LineWidth',1.4); hold on;
    plot(cmd.Time, cmd.Data, 'r--','LineWidth',1.0);
    grid on; box on;
    title(sprintf('\\omega = %.1f rad/sec',w),'FontWeight','bold');
    ylabel('pitch rate  [deg/s]');
    legend({'q  (output)','K\cdotq_c  (command)'},'Location','best','FontSize',8);
    if w > 5.1, ylim([-400 400]); end

    % ---- bottom row: rate limiter element ----
    nexttile(j+2);
    plot(dc.Time, dc.Data, 'm-','LineWidth',1.0); hold on;
    plot(dlt.Time,dlt.Data,'k-','LineWidth',1.4);
    grid on; box on;
    xlabel('Time (sec)'); ylabel('RLE  [deg]');
    legend({'\delta_c  (commanded)','\delta  (rate limited)'}, ...
           'Location','best','FontSize',8);
    if w > 5.1, ylim([-400 400]); end
end

title(tl,'Figure 2-8 (Simulink reproduction):  Time Domain Analysis', ...
      'FontWeight','bold','FontSize',12);

fprintf('Model "%s" built and simulated. Leave it open to inspect/edit.\n',mdl);
