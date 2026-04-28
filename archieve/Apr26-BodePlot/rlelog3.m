%% Rate-Limiter Describing Function via Sine Sweep
% ------------------------------------------------------------------------
%  Drives a Simulink rate-limiter element with a continuous sinusoid at
%  each frequency in the sweep and extracts the first-harmonic describing
%  function:
%
%      N(A,w) = ( b1 + j*a1 ) / A
%
%  where, for input r(t) = A*sin(w*t) and steady-state output y(t):
%      b1 = (2/T) * Integral( y(t) * sin(w*t) dt )   over the last N_capture periods
%      a1 = (2/T) * Integral( y(t) * cos(w*t) dt )   over the last N_capture periods
%
%  Model:  rle2.slx
%      Sine Wave (Amplitude=A, Frequency=w)  ->  Rate Limiter (+/-R)  ->  yout
%
%  Changes vs. rlelog2.m:
%    1. G is preallocated (was growing inside the loop).
%    2. Stale top-level "t_stop = 30" removed.
%    3. Transient is dropped EXPLICITLY: integrate over the last N_capture
%       periods (default 3) of N_total simulated periods (default 10),
%       not just the very last period.
%    4. Trivial regime w <= R/A is short-circuited (N = 1 + 0j) - saves
%       wall time on the low-frequency end where the rate limiter is
%       provably inactive.
%    5. load_system used so the model is already compiled before the loop.
% ------------------------------------------------------------------------

clear; clc; close all
tic

%% Generating Sinwe Wave

%
t_stop = 30;      % simulation stop time
num_t_step = 1000;      % numbers of time step
t = linspace(0, t_stop, num_t_step)';      % time
X_i_max = 1;
A = X_i_max;
%
R = 0.725;        % slew rate
%
modelName = 'rle'; % Rate Limiter Element

u = zeros(num_t_step,10);
cnt = 0;

% omega = 0.1:0.1:10; % Linear increments
omega = logspace(-1,2,100); % logarithmic increments

for i = 1:length(omega)
    w = omega(i);
    u(:, i) = X_i_max*sin(w*t);        % input signal
end

inputData = [t, u];          % [time value] format for From Workspace
assignin('base','myInput',inputData);
% In the model, set From Workspace block Data = myInput
out = sim(modelName);
% After sim finishes, read output variable from To Workspace, e.g.:
outData = out.yout;           % if To Workspace used default name
t = out.tout;



%% Parameters
A         = 1;                        % input amplitude
R         = 0.725;                    % rate limit  (same units as A per second)
freqs     = logspace(-1, 2, 100);     % sweep frequencies (rad/s)
N_total   = 10;                       % periods simulated per frequency
N_capture = 3;                        % final periods used for Fourier integration
modelName = 'rle3';

assert(N_capture <= N_total, 'N_capture must not exceed N_total.');

%% Pre-allocate
G = zeros(size(freqs));               % complex describing function

%% Sweep
load_system(modelName);
w_onset = R / A;                      % rate-saturation onset frequency

for i = 1:length(freqs)
    w        = freqs(i);              % read by Sine Wave block at sim start
    T_period = 2*pi/w;
    t_stop   = N_total * T_period;

    % Below saturation onset the rate limiter cannot activate (peak input
    % rate A*w <= R), so y == r and N = 1.  Skip the sim.
    if w <= w_onset
        G(i) = 1 + 0i;
        continue
    end

    simOut = sim(modelName, 'StopTime', num2str(t_stop));
    t = simOut.tout;
    y = simOut.yout(:,1);

    % Integrate over the LAST N_capture periods (transient discarded)
    t_window = N_capture * T_period;
    idx      = t >= (t_stop - t_window);
    tc       = t(idx);
    yc       = y(idx);

    % First-harmonic Fourier coefficients
    Tc = tc(end) - tc(1);
    b1 = (2/Tc) * trapz(tc, yc .* sin(w*tc));
    a1 = (2/Tc) * trapz(tc, yc .* cos(w*tc));

    G(i) = (b1 + 1i*a1) / A;
end

mag_dB = 20*log10(abs(G));
phase  = rad2deg(angle(G));

%% FRD object (handy for series interconnections downstream)
sys = frd(G, freqs);

%% Plots
figure('Name','Rate-Limiter Describing Function','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
semilogx(freqs, mag_dB, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.4); grid on
ylabel('Magnitude  (dB)');
title('Rate-Limiter Describing Function — Magnitude');
xline(w_onset, 'r--', 'Label','\omega = R/A', 'LabelOrientation','horizontal');

ax2 = nexttile;
semilogx(freqs, phase, 'r-o', 'MarkerSize', 4, 'LineWidth', 1.4); grid on
ylabel('Phase  (deg)'); xlabel('\omega  (rad/s)');
title('Rate-Limiter Describing Function — Phase');
xline(w_onset, 'r--', 'Label','\omega = R/A', 'LabelOrientation','horizontal');
yline(-90, 'k--', 'Label','-90\circ', 'LabelOrientation','horizontal');

linkaxes([ax1 ax2],'x');
xlim([min(freqs) max(freqs)]);

toc
