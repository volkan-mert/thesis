%% mehra_fig7_fig8.m
% Reproduce Figures 7 & 8 of:
%   Mehra, R., & Prasanth, R. (1998). "Bifurcation and limit cycle analysis
%   of nonlinear pilot induced oscillations." AIAA-98-4249.
%
% Model: Mehra.slx
%   Sin(qco, w) -> Sum(+-) -> Pilot Kp/(thau*s+1) -> RateLimiter(+/-R)
%   -> Airframe [86.9 79.2528 2.2252]/[1 2.06 5.94 2.03 0.053] -> theta (fb)
%
% NOTE: the pilot lag variable in the model is named THAU (not tau).
% NOTE: block "Airframe theta/delta" contains '/', escaped as '//' in paths.

clear; clc; close all

%% ---------------- User settings ----------------
mdl         = "Mehra";
Tstop       = 80;            % [s] long enough to settle onto the limit cycle
ssFrac      = 0.70;          % amplitude measured over the last 30% of the run
useParallel = true;          % parsim if Parallel Computing Toolbox is present

% System parameters
R    = 15;                   % slew-rate limit [deg/s]
qco  = 1;                    % pilot command (kick) amplitude

% Sweep grids
fHzList = 0.1:0.1:2;                          % command band [0.1 2] Hz
wList   = 2*pi*fHzList;                       % [rad/s]
tauList = [0 0.03 0.06 0.09];                 % four curves of Fig. 7
KpList  = [0.1:0.1:2.4, 2.45:0.05:3.4, ...   % refined near the jump (~2.7)
           3.6:0.2:5];

%% ---------------- One-time model preparation ----------------
load_system(mdl);

% Enable signal logging on the airframe output (theta)
phAir = get_param(mdl + "/Airframe theta//delta", "PortHandles");
set_param(phAir.Outport(1), "DataLogging","on", ...
    "DataLoggingNameMode","Custom", "DataLoggingName","theta");

set_param(mdl, "StopTime", num2str(Tstop), ...
               "SignalLogging","on", "SignalLoggingName","logsout", ...
               "SaveOutput","off", "SaveState","off");
save_system(mdl);            % logging config must be saved for Fast Restart

%% ---------------- Build the SimulationInput array ----------------
[TAU, KP, W] = ndgrid(tauList, KpList, wList);   % size: [nTau x nKp x nW]
n = numel(TAU);
in(1:n) = Simulink.SimulationInput(mdl);
for i = 1:n
    in(i) = in(i) ...
        .setVariable("thau", TAU(i)) ...   % <-- model variable is 'thau'
        .setVariable("Kp",   KP(i)) ...
        .setVariable("w",    W(i)) ...
        .setVariable("qco",  qco) ...
        .setVariable("R",    R);
end

%% ---------------- Run the sweep ----------------
fprintf("Running %d simulations...\n", n);
if useParallel && license("test","Distrib_Computing_Toolbox")
    out = parsim(in, "ShowProgress","on", "UseFastRestart","on", ...
                     "TransferBaseWorkspaceVariables","off");
else
    out = sim(in, "ShowProgress","on", "UseFastRestart","on");
end

%% ---------------- Extract steady-state limit-cycle amplitude ----------------
Amp = nan(size(TAU));
for i = 1:n
    if ~isempty(out(i).ErrorMessage)
        warning("Run %d failed: %s", i, out(i).ErrorMessage);
        continue
    end
    ts  = out(i).logsout.getElement("theta").Values;
    t   = ts.Time;  y = squeeze(ts.Data);
    idx = t >= ssFrac*t(end);
    Amp(i) = (max(y(idx)) - min(y(idx)))/2;   % half peak-to-peak
end

% Fig. 7 shows the limit-cycle amplitude over the [0.1 2] Hz band:
% take the worst case (max) across the command-frequency dimension.
AmpKp = max(Amp, [], 3);                      % [nTau x nKp]

%% ---------------- Figure 7: amplitude vs pilot gain ----------------
figure("Name","Mehra Fig. 7","Color","w"); hold on; box on
styles = {'--','-.',':','-'};
for m = 1:numel(tauList)
    plot(KpList, AmpKp(m,:), styles{m}, "Color","k", "LineWidth",1.1, ...
        "DisplayName", sprintf("\\tau = %.2f", tauList(m)));
end
xlabel("Pilot Gain"); ylabel("Limit Cycle Amplitude");
title({"X-15 limit cycle amplitude vs. pilot gain", ...
       "(jump phenomenon; PIO for K_p > \approx 2.7)"});
legend("Location","northwest"); xlim([0 5]);

%% ---------------- Figure 8: phase portrait at Kp = 5 ----------------
% Single long run; theta_dot computed numerically (the Integrator block
% feeding the XY Graph actually outputs \int(theta), not d(theta)/dt).
in8 = Simulink.SimulationInput(mdl) ...
      .setVariable("thau", 0.09) ...
      .setVariable("Kp",   5) ...
      .setVariable("w",    2*pi*0.5) ...
      .setVariable("qco",  qco) ...
      .setVariable("R",    R) ...
      .setModelParameter("StopTime","120");
out8 = sim(in8);

ts8 = out8.logsout.getElement("theta").Values;
t8  = ts8.Time;  th8 = squeeze(ts8.Data);
thd8 = gradient(th8, t8);                     % numerical d(theta)/dt

figure("Name","Mehra Fig. 8","Color","w");
plot(th8, thd8, "k", "LineWidth", 0.5); grid on; box on
xlabel("\theta"); ylabel("\theta_{dot}");
title("X-15 Limit Cycle at K_p = 5");
axis tight

%% ---------------- Save results ----------------
save("mehra_sweep_results.mat", "Amp","AmpKp","TAU","KP","W", ...
     "tauList","KpList","fHzList","R","qco","Tstop","ssFrac");
fprintf("Done. Results saved to mehra_sweep_results.mat\n");
