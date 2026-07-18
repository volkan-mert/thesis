%% draw_figure7.m
%  Reproduce Mehra "Figure 7": X-15 limit-cycle amplitude vs pilot gain Kp,
%  for tau = {0, 0.03, 0.06, 0.09}, by driving the Simulink model Mehra.slx.
%
%  WHAT THIS DOES
%    - Loops over (tau, Kp), sets the base-workspace variables the blocks read
%      (Kp, thau, R, qco, f), runs the model, and measures the steady-state
%      peak amplitude of theta -> the (time-domain) limit-cycle amplitude.
%
%  TWO IMPORTANT NOTES (read the response that came with this file):
%    1) ACTUATOR. In the delivered model the actuator is wired as
%          Sum(u - delta) -> [built-in Rate Limiter] -> delta
%       In the unsaturated region that resolves to a STATIC gain of 0.5 with
%       NO phase lag (and it forms an algebraic loop). That is not the Fig.6
%       actuator, and with it the tau=0 loop never destabilises -> no PIO.
%       Set DO_FIX_ACTUATOR = true below to build a corrected copy that uses
%       the Fig.6 actuator:  delta_dot = sat( 20*(u - delta), +/-R ), i.e.
%       first-order 20/(s+20) whose RATE is limited to +/-15 deg/s.
%    2) THE JUMP. Fig.7's discontinuous "jump" is a harmonic-balance / BACTM
%       (describing-function continuation) result. A forward time simulation
%       settles onto a single limit-cycle branch and gives a smooth curve, not
%       the printed jump. Use this script to get the time-domain amplitude
%       curves; use the describing-function method (your gilbreathDF/OLOP
%       tooling) to reproduce the exact jump.
%
%  Author helper: run section by section or all at once.

clear; clc;

%% ---------------- configuration ----------------
srcModel        = 'Mehra';          % base model (Mehra.slx must be on the path)
DO_FIX_ACTUATOR = true;             % true -> build & simulate the Fig.6 actuator
USE_PARSIM      = false;            % true -> Parallel Computing Toolbox (parsim)

R    = 15;        % slew-rate limit [deg/s]
qco  = 1;         % pilot command (kick) amplitude
fHz  = 0.5;       % command frequency [Hz], within the paper's [0.1 2] band

tauList = [0 0.03 0.06 0.09];       % the four curves
KpList  = 0.1:0.1:5;                % pilot-gain grid (refine near the jump)

Tstop   = 300;    % [s] long enough for the slow airframe pole (~35 s time const)
measWin = 60;     % [s] final window used to measure the limit-cycle amplitude

%% ---------------- prepare the model ----------------
load_system(srcModel);

if DO_FIX_ACTUATOR
    model = buildFixedActuator(srcModel);   % returns 'Mehra_fixed'
else
    model = srcModel;
end

% Log the airframe output signal "theta".  NOTE the block name contains a
% literal '/', which is escaped as '//' inside a Simulink path.
afBlock = [model '/Airframe theta//delta'];
ph = get_param(afBlock,'PortHandles');
set_param(ph.Outport(1),'DataLogging','on', ...
                        'DataLoggingNameMode','Custom', ...
                        'DataLoggingName','theta');
set_param(model,'SignalLogging','on','SignalLoggingName','logsout');

%% ---------------- build the simulation input set ----------------
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
        in = in.setVariable('qco',  qco);
        in = in.setVariable('f',    fHz);
        in = in.setModelParameter('StopTime', num2str(Tstop));
        in = in.setModelParameter('SignalLogging','on');
        % A stiff/implicit solver handles the near-integrator + saturation well:
        in = in.setModelParameter('Solver','ode23t', ...
                                  'RelTol','1e-6','MaxStep','0.02');
        inList(idx) = in;
    end
end

%% ---------------- run ----------------
if USE_PARSIM
    outList = parsim(inList,'ShowProgress','on');
else
    outList = sim(inList);
end

%% ---------------- measure amplitude ----------------
amp = zeros(nTau,nKp);
idx = 0;
for i = 1:nTau
    for j = 1:nKp
        idx = idx + 1;
        ts = getThetaSignal(outList(idx));
        t  = ts.Time;  y = ts.Data;
        seg = y(t >= (Tstop - measWin));
        amp(i,j) = 0.5*(max(seg) - min(seg));   % 0.5*peak-to-peak
    end
end

%% ---------------- plot (Figure 7 style) ----------------
figure('Color','w'); hold on; grid on; box on;
sty = {'--','-.',':','-'};
for i = 1:nTau
    plot(KpList, amp(i,:), sty{i}, 'LineWidth', 1.6, 'Color','k');
end
xlabel('Pilot Gain'); ylabel('Limit Cycle Amplitude');
legend(compose('\\tau = %.2f', tauList), 'Location','northwest');
title('X-15 Limit Cycle Amplitude vs Pilot Gain');
xlim([0 5]);

%% ================= helper functions =================
function ts = getThetaSignal(out)
% Robustly pull the logged "theta" signal from a SimulationOutput.
    ls = out.logsout;
    try
        ts = ls.get('theta').Values;
    catch
        ts = ls{1}.Values;    % fallback: first logged signal
    end
end

function newModel = buildFixedActuator(srcModel)
% Create Mehra_fixed with the Fig.6 actuator:
%   delta_dot = sat( 20*(u - delta), +/-R )   (20/(s+20) with rate limit R)
% This replaces the built-in "Rate Limiter" block. The Sum (u - delta) and
% the delta feedback already exist in the model; we insert Gain(20) ->
% Saturation([-R R]) -> Integrator(1/s), and rewire delta.
    newModel = [srcModel '_fixed'];
    if bdIsLoaded(newModel), close_system(newModel,0); end
    save_system(srcModel, newModel, 'OverwriteIfChanged', true);
    load_system(newModel);

    rl  = [newModel '/Rate Limiter'];
    sumB= [newModel '/Sum'];                    % computes delta_c = u - delta
    af  = [newModel '/Airframe theta//delta'];

    % Remember Rate Limiter position, then delete it and the lines touching it.
    pos = get_param(rl,'Position');
    delete_line_ifexists(newModel, 'Sum/1', 'Rate Limiter/1');
    % delta feedback (Rate Limiter -> Sum in2) and forward (Rate Limiter -> Airframe)
    delete_line_ifexists(newModel, 'Rate Limiter/1', 'Sum/2');
    delete_line_ifexists(newModel, 'Rate Limiter/1', 'Airframe theta//delta/1');
    delete_block(rl);

    % New actuator blocks (placed left-to-right where the Rate Limiter was).
    g = [newModel '/ActGain'];  s = [newModel '/RateSat'];  ii = [newModel '/ActInt'];
    add_block('simulink/Math Operations/Gain', g, ...
        'Gain','20','Position',[pos(1)      pos(2)   pos(1)+30 pos(4)]);
    add_block('simulink/Discontinuities/Saturation', s, ...
        'UpperLimit','R','LowerLimit','-R', ...
        'Position',[pos(1)+55  pos(2)   pos(1)+85 pos(4)]);
    add_block('simulink/Continuous/Integrator', ii, ...
        'Position',[pos(1)+110 pos(2)   pos(1)+140 pos(4)]);

    % Wire:  Sum -> Gain -> Sat -> Integrator(=delta) -> {Airframe, Sum in2}
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
