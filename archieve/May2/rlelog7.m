clear;clc;close all % clean start
%
tic
%
t_fss = 1e-3; % fized step size
%
%% Example : Draw Fig. 16 Nichols chart: open-loop frequency response q∕qc 
% 2025, Rodrigues et al, CEAS Aeronautical Journal (2026), 
% "Conditional integrator sliding mode control to reduce susceptibility 
% to pilot‑induced oscillations"

num_claw_7 = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7 = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7 = tf(num_claw_7, den_claw_7);

num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7 = tf(num_ldynac_7, den_ldynac_7);

Gs_ac_7 = Gs_claw_7*Gs_ldynac_7; % TF of the open-loop by q / q_c.

Gs_ac = Gs_ac_7;

[num_ac, den_ac] = tfdata(Gs_ac, 'v'); % This is unstable system that is not a Bounded Input and Bounded Output System

Kp = 1; % Pilot Gain

pick_cursor_magnitude_olop = 2.0629; % omega_OLOP Magnitude in dB
pick_cursor_phase_olop = -136.6216; % omega_OLOP in degrees

figure(99);
nichols(Gs_ac);
opts = nicholsoptions;
opts.PhaseMatching      = 'on';
opts.PhaseMatchingFreq  = 1;        % rad/s
opts.PhaseMatchingValue = -180;     % anchor phase at that frequency
nicholsplot(Gs_ac, opts); 
hold on
xline(-180, '--'); % in degrees0
yline(0, '--'); % in dBs
hold on
plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop,'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
hold on
text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, '\omega_{OLOP} = 5.1', 'Color', 'red', 'FontWeight', 'bold', 'FontSize', 12);
%
grid on
xlim([-300 -50]);
ylim([-20 15]);

title('Nichols Chart of Aircraft and Control Law Transfer Functions Only');

%% The Open-loop Onset Point Method

%
num_t_step = 1e4;      % numbers of time step (per frequency)
X_i_max = 1;
A = X_i_max;
%
R = 0.725;        % slew rate
%
modelName1 = 'rle7sinewave'; % sine wave generation
modelName2 = 'rle7';         % Rate Limiter + Aircraft + Pilot loop

load_system(modelName1);
load_system(modelName2);

% Open scopes (note: scope will only show the LAST sweep iteration)
scopes = find_system(modelName2, 'BlockType', 'Scope');
for sc = 1:length(scopes)
    open_system(scopes{sc});
    cfg = get_param(modelName2 + "/Scope", 'ScopeConfiguration');

    % Fix for the Y-axis autoscaling error
    cfg.AxesScaling = 'Auto';

    % Proactive fix for the X-axis (time) autoscaling
    cfg.TimeSpan = 'Auto';

    % Add this line to show the legend
    cfg.ShowLegend = true;
end

%% Frequency sweep — scalar-input loop (fixes the SISO Aircraft Block dim error)

omega = logspace(-1, 2, 100);     % logarithmic increments
H     = zeros(1, length(omega));  % preallocate first-harmonic frequency response

for i = 1:length(omega)
    w        = omega(i);
    T_period = 2*pi/w;
    t_stop   = 10 * T_period;                       % 10 full periods at this omega
    t        = linspace(0, t_stop, num_t_step)';

    % ---- 1) Generate sinusoidal input via rle6sinewave -----------------
    %      (rle6sinewave reads X_i_max and w from base workspace)
    out1 = sim(modelName1, 'StopTime', num2str(t_stop));
    u_i  = interp1(out1.myinput.Time, out1.myinput.Data, t, 'linear', 'extrap');

        if any(isnan(u_i)) || any(isinf(u_i))
            error('Input signal contains Inf/NaN at omega = %.4f rad/s!', w);
        % else
        %     disp('Input signal doesn''t contain any Inf/NaN. It''s all right.');
        end

    % ---- 2) Push SCALAR-channel timeseries to base workspace -----------
    %      The From Workspace block in rle6 reads `myInput`.
    inputDataTS = timeseries(u_i, t);
    assignin('base', 'myInput', inputDataTS);

    % ---- 3) Run the closed loop (Rate Limiter -> Aircraft -> Pilot) ----
    out2    = sim(modelName2, 'StopTime', num2str(t_stop));
    y       = out2.yout;
    t_out   = out2.tout;

    % ---- 4) First-harmonic Fourier coefficients over the last period ---
    idx     = t_out > (t_out(end) - T_period);
    t_last  = t_out(idx);
    y_last  = y(idx);

    b1 = (2/T_period) * trapz(t_last, y_last .* sin(w * t_last));
    a1 = (2/T_period) * trapz(t_last, y_last .* cos(w * t_last));

    H(i) = (b1 + 1i*a1) / A;

    if mod(i,10) == 0
        fprintf('  swept %3d / %3d  (omega = %8.4f rad/s)\n', ...
                 i, length(omega), w);
    end
   
end

mag   = abs(H);
phase = rad2deg(angle(H));

%% Frequency Response Data-model

sys_frd = frd(H, omega);

% sys_est = tfest(sys_frd, 1, 0);
% 
% [num_ac, den_ac] = tfdata(sys_est, 'v');

%% Nichols Chart via FRD
hold on
f6 = figure('Name', 'Nichols Chart (First Harmonic) by FRD', 'NumberTitle', 'off');
nichols(sys_frd);
grid on;

toc
