%% draw_figure7_chirp.m
%  Draw the "Figure 7" family (limit-cycle amplitude vs pilot gain Kp, for
%  tau = {0,0.03,0.06,0.09}) by driving Mehra.slx with a CHIRP command that
%  sweeps 0.1 -> 2 Hz, then sweeping Kp and measuring the peak theta response.
%
%  It (optionally) builds a corrected copy of the model with:
%    - the Fig.6 actuator  delta_dot = sat( 20*(u-delta), +/-R )   [20/(s+20)]
%    - the sine "Nonlinear Oscillator" replaced by a Chirp Signal block
%
%  READ THIS FIRST -----------------------------------------------------------
%  This loop is LINEARLY UNSTABLE for Kp > ~0.3 (crossing at ~5.4 rad/s), so a
%  rate-limited limit cycle exists across the whole Kp range. A forward-time
%  chirp sweep therefore produces a SMOOTH, tau-ordered curve -- NOT the
%  discontinuous jump at Kp~2.7 in the printed figure. That jump is a fold in
%  the harmonic-balance (describing-function / BACTM) branch and can only be
%  reproduced by that continuation method, not by simulation. Use this script
%  as a broadband PIO-excitation diagnostic; use the DF method for the exact
%  Fig.7.
%  ---------------------------------------------------------------------------

clear; clc;

%% ---------------- configuration ----------------
srcModel        = 'Mehra';
DO_FIX_ACTUATOR = true;      % Fig.6 actuator (needed or tau=0 never limit-cycles)
USE_PARSIM      = false;

R       = 15;                % slew-rate limit [deg/s]
f_start = 0.1;              % chirp initial frequency [Hz]
f_end   = 2.0;              % chirp final frequency  [Hz]

tauList = [0 0.03 0.06 0.09];
KpList  = 0.1:0.1:5;

Tsweep  = 120;              % [s] chirp duration = StopTime (slow sweep -> good dwell)
t_skip  = 5;               % [s] ignore startup transient before measuring

%% ---------------- prepare the model ----------------
load_system(srcModel);
model = srcModel;
if DO_FIX_ACTUATOR, model = buildFixedActuator(model); end
model = swapInChirp(model, f_start, f_end, Tsweep);

% Log airframe output "theta" (block name has a literal '/', escaped as '//').
afBlock = [model '/Airframe theta//delta'];
ph = get_param(afBlock,'PortHandles');
set_param(ph.Outport(1),'DataLogging','on', ...
                        'DataLoggingNameMode','Custom','DataLoggingName','theta');
set_param(model,'SignalLogging','on','SignalLoggingName','logsout');

%% ---------------- build simulation inputs ----------------
nTau = numel(tauList);  nKp = numel(KpList);
inList(1,nTau*nKp) = Simulink.SimulationInput(model);
idx = 0;
for i = 1:nTau
    for j = 1:nKp
        idx = idx + 1;
        in = Simulink.SimulationInput(model);
        in = in.setVariable('Kp',   KpList(j));
        in = in.setVariable('thau', tauList(i));
        in = in.setVariable('R',    R);
        in = in.setVariable('qco',  1);     % harmless if unused by chirp
        in = in.setVariable('f',    f_start);
        in = in.setModelParameter('StopTime', num2str(Tsweep));
        in = in.setModelParameter('SignalLogging','on');
        in = in.setModelParameter('Solver','ode23t','RelTol','1e-6','MaxStep','0.02');
        inList(idx) = in;
    end
end

%% ---------------- run ----------------
if USE_PARSIM
    outList = parsim(inList,'ShowProgress','on');
else
    outList = sim(inList);
end

%% ---------------- measure peak amplitude over the sweep ----------------
amp = zeros(nTau,nKp);
idx = 0;
for i = 1:nTau
    for j = 1:nKp
        idx = idx + 1;
        ts = getTheta(outList(idx));
        t = ts.Time; y = ts.Data;
        seg = y(t >= t_skip);
        amp(i,j) = 0.5*(max(seg) - min(seg));   % half peak-to-peak
    end
end

%% ---------------- plot ----------------
figure('Color','w'); hold on; grid on; box on;
sty = {'--','-.',':','-'};
for i = 1:nTau
    plot(KpList, amp(i,:), sty{i}, 'LineWidth', 1.6, 'Color','k');
end
xlabel('Pilot Gain'); ylabel('Limit Cycle Amplitude (peak \theta)');
legend(compose('\\tau = %.2f', tauList),'Location','northwest');
title('Chirp-excited (0.1\rightarrow2 Hz) amplitude vs pilot gain');
xlim([0 5]);

%% ================= helpers =================
function ts = getTheta(out)
    ls = out.logsout;
    try, ts = ls.get('theta').Values; catch, ts = ls{1}.Values; end
end

function m = swapInChirp(model, f1, f2, Tsweep)
% Replace the sine "Nonlinear Oscillator" with a Chirp Signal block that
% sweeps f1 -> f2 over [0, Tsweep]. (Chirp Signal output amplitude is 1;
% add a Gain block after it if you need a different command amplitude.)
    m = model;
    sine = [m '/Nonlinear Oscillator'];
    pos  = get_param(sine,'Position');
    try, delete_line(m,'Nonlinear Oscillator/1','Sum1/1'); catch, end
    delete_block(sine);
    ch = [m '/Chirp'];
    add_block('simulink/Sources/Chirp Signal', ch, 'Position', pos, ...
        'f1', num2str(f1), 'T', num2str(Tsweep), 'f2', num2str(f2));
    add_line(m,'Chirp/1','Sum1/1','autorouting','on');
    save_system(m);
end

function newModel = buildFixedActuator(srcModel)
% Fig.6 actuator: delta_dot = sat( 20*(u-delta), +/-R ) replacing the
% built-in Rate Limiter. The Sum (u-delta) and delta feedback already exist.
    newModel = [srcModel '_fixed'];
    if bdIsLoaded(newModel), close_system(newModel,0); end
    save_system(srcModel, newModel, 'OverwriteIfChanged', true);
    load_system(newModel);

    rl = [newModel '/Rate Limiter'];
    pos = get_param(rl,'Position');
    delete_line_ifexists(newModel,'Sum/1','Rate Limiter/1');
    delete_line_ifexists(newModel,'Rate Limiter/1','Sum/2');
    delete_line_ifexists(newModel,'Rate Limiter/1','Airframe theta//delta/1');
    delete_block(rl);

    g  = [newModel '/ActGain'];  s = [newModel '/RateSat'];  ii = [newModel '/ActInt'];
    add_block('simulink/Math Operations/Gain', g, ...
        'Gain','20','Position',[pos(1)      pos(2) pos(1)+30  pos(4)]);
    add_block('simulink/Discontinuities/Saturation', s, ...
        'UpperLimit','R','LowerLimit','-R', ...
        'Position',[pos(1)+55  pos(2) pos(1)+85  pos(4)]);
    add_block('simulink/Continuous/Integrator', ii, ...
        'Position',[pos(1)+110 pos(2) pos(1)+140 pos(4)]);

    add_line(newModel,'Sum/1','ActGain/1','autorouting','on');
    add_line(newModel,'ActGain/1','RateSat/1','autorouting','on');
    add_line(newModel,'RateSat/1','ActInt/1','autorouting','on');
    add_line(newModel,'ActInt/1','Airframe theta//delta/1','autorouting','on');
    add_line(newModel,'ActInt/1','Sum/2','autorouting','on');
    save_system(newModel);
end

function delete_line_ifexists(mdl, src, dst)
    try, delete_line(mdl, src, dst); catch, end
end
