clear;clc;close all % clean start
%
tic
%% Example 1: Aircraft Block — X-15 PIO during flare

% (Amato, Iervolino, Pandit, Scala, Verde, CDC 2000: "Analysis of Pilot-in-the-Loop Oscillations Due to Position and Rate Saturations")
num_ac_1 = [3.476 3.1708072 0.0896237936];
den_ac_1 = [1 0.8608 5.3159942 0.108928 0.0529];
Gs_ac_1 = tf(num_ac_1, den_ac_1);

% Pilot Gain (20 deg PM of the linearised model)
Kp_ac_1 = 3.5;

%% Example 2: Flight control system design considering rate saturation 

% (1987, Duda, DLR, "Flight control system design considering rate saturation" for Kp = 0.68)
num_ac_2 = [4066.3, 17509.4868];
den_ac_2 = [1, 44.85, 1023.0574, 786.6747, -2038.062];
Gs_ac_2 = tf(num_ac_2, den_ac_2);

Kp_ac_2 = 0.68; % OPEN-LOOP TF — Low-gain design (K_fb = 0.68)

%% Example 3: Flight control system design considering rate saturation 

% (1987, Duda, DLR, "Flight control system design considering rate saturation" for Kp = 1.0)
num_ac_3 = [5979.8, 25748.8188];
den_ac_3 = [1, 44.85, 1023.0574, 786.6747, -2038.062];
Gs_ac_3 = tf(num_ac_3, den_ac_3);

Kp_ac_3 = 1; % OPEN-LOOP TF — High-gain design (K_fb = 0.68)

%% Example 4: Pilot-Induced Oscillations and Their Prevention

% (2020, Adrievsky, CPHS, "Pilot-Induced Oscillations and Their Prevention" for Kp = 2.3)
num_ac_4 = [86.9, 79.27, 2.24];
den_ac_4 = [1, 27.06, 57.45, 150.75, 50.73, 1.32];
Gs_ac_4 = tf(num_ac_4, den_ac_4);

Kp_ac_4 = 2.3; % K_fb = 2.3

%% Example 5: Pilot-Induced Oscillations and Their Prevention

% (2020, Adrievsky, CPHS, "Pilot-Induced Oscillations and Their Prevention" for Kp = 2.7)
num_ac_5 = [86.9, 79.27, 2.24];
den_ac_5 = [1, 27.06, 57.45, 150.75, 50.73, 1.32];
Gs_ac_5 = tf(num_ac_5, den_ac_5);

Kp_ac_5 = 2.7; % K_fb = 2.7

%% Example 6: Simple Pitch-Attitude Autopilot of F-16
% % (2010, F. Lewis, Book, Aircraft Control and Simulation, Chapter 4: Autopilots, Example 4.6: Simple Pitch-Attitude Autopilot)

num_ac_6 = [-1.133, -0.65242, -0.01207, -1.9934e-06];
den_ac_6 = [1.0, 1.05193, 1.76849, 0.01747, 0.01419, 2.6845e-06];
Gs_ac_6 = tf(num_ac_6, den_ac_6);

Kp_ac_6 = 4; % K_fb = 4

%% Deployment of the Parameters of the Blocks of Simulink Model

% % TF of the paper of Amato et al 2000, CDC
% num_ac = num_ac_1; % numerator of the transfer function of the aircraft
% den_ac = den_ac_1; % denominator of the transfer function of the aircraft
% Gs_ac = Gs_ac_1;  % transfer function of the aircraft
% Kp = Kp_ac_1; % pilot of the transfer function of the aircraft 

% TF of Duda, 1987, Flight control system design considering rate saturation
% num_ac = num_ac_2; % numerator of the transfer function of the aircraft
% den_ac = den_ac_2; % denominator of the transfer function of the aircraft
% Gs_ac = Gs_ac_2;  % transfer function of the aircraft
% Kp = Kp_ac_2; % pilot of the transfer function of the aircraft 

% (1987, Duda, DLR, "Flight control system design considering rate saturation" for Kp = 1.0)
% num_ac = num_ac_3; % numerator of the transfer function of the aircraft
% den_ac = den_ac_3; % denominator of the transfer function of the aircraft
% Gs_ac = Gs_ac_3;  % transfer function of the aircraft
% Kp = Kp_ac_3; % pilot of the transfer function of the aircraft

% (2020, Adrievsky, CPHS, "Pilot-Induced Oscillations and Their Prevention" for Kp = 2.3)
% num_ac = num_ac_4; % numerator of the transfer function of the aircraft
% den_ac = den_ac_4; % denominator of the transfer function of the aircraft
% Gs_ac = Gs_ac_4;  % transfer function of the aircraft
% Kp = Kp_ac_4; % pilot of the transfer function of the aircraft

% % (2020, Adrievsky, CPHS, "Pilot-Induced Oscillations and Their Prevention" for Kp = 2.7)
% num_ac = num_ac_5; % numerator of the transfer function of the aircraft
% den_ac = den_ac_5; % denominator of the transfer function of the aircraft
% Gs_ac = Gs_ac_5;  % transfer function of the aircraft
% Kp = Kp_ac_5; % pilot of the transfer function of the aircraft

% % (2010, F. Lewis, Book, Aircraft Control and Simulation, Chapter 4: Autopilots, Example 4.6: Simple Pitch-Attitude Autopilot)
% num_ac = num_ac_6; % numerator of the transfer function of the aircraft
% den_ac = den_ac_6; % denominator of the transfer function of the aircraft
% Gs_ac = Gs_ac_6;  % transfer function of the aircraft
% Kp = Kp_ac_6; % pilot of the transfer function of the aircraft

f0 = figure('Name', 'Bode Plot of Aircraft Transfer Function only', 'NumberTitle', 'off');
bode(Kp*Gs_ac);
grid on
title('Bode Plot of Aircraft Transfer Function only');

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
% scopes = find_system(modelName2, 'BlockType', 'Scope');
% for sc = 1:length(scopes)
%     open_system(scopes{sc});
%     cfg = get_param(modelName2 + "/Scope", 'ScopeConfiguration');
% 
%     % Fix for the Y-axis autoscaling error
%     cfg.AxesScaling = 'Auto';
% 
%     % Proactive fix for the X-axis (time) autoscaling
%     cfg.TimeSpan = 'Auto';
% 
%     % Add this line to show the legend
%     cfg.ShowLegend = true;
% end

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

%% Bode Plot (manual)
f1 = figure('Name', 'Nonlinear Bode Plot (First Harmonic)', 'NumberTitle', 'off');
subplot(2,1,1);
semilogx(omega, 20*log10(mag), 'b-o', 'LineWidth', 1.5);
grid on;
ylabel('Magnitude (dB)');
title('Nonlinear Bode Plot (First Harmonic)');
xline(R/A,'r--','Label','R/A','LabelOrientation','horizontal');

subplot(2,1,2);
semilogx(omega, phase, 'r-o', 'LineWidth', 1.5);
grid on;
ylabel('Phase (deg)');
xlabel('Frequency (rad/s)');
xline(R/A,'r--','Label','R/A','LabelOrientation','horizontal');

%% Bode Plot via FRD
f2 = figure('Name', 'Nonlinear Bode Plot (First Harmonic) by FRD', 'NumberTitle', 'off');
bode(sys_frd);
hold on;
xline(R/A,'r--','Label','R/A','LabelOrientation','horizontal');
grid on;

%% Nyquist (manual)
f3 = figure('Name', 'Nyquist Diagram (First Harmonic)', 'NumberTitle', 'off');
plot(real(H),  imag(H), 'b',  'LineWidth', 1.5); hold on;
plot(real(H), -imag(H), 'b--','LineWidth', 1.5);   % mirror for -ω
plot(-1, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal;
xlabel('Real Axis'); ylabel('Imaginary Axis');
title('Nyquist Diagram');

%% Nyquist via FRD
f4 = figure('Name', 'Nyquist Diagram (First Harmonic) by FRD', 'NumberTitle', 'off');
nyquist(sys_frd);
grid on;

%% Nichols Chart (manual)
f5 = figure('Name', 'Nichols Chart (First Harmonic)', 'NumberTitle', 'off');
plot(phase, 20*log10(mag), 'b', 'LineWidth', 1.5); hold on;   % dB on y-axis
plot(-180, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Open-Loop Phase (deg)');
ylabel('Open-Loop Gain (dB)');
title('Nichols Chart');
grid on;
% xlim([-270 -90]); ylim([-40 20]);

%% Nichols Chart via FRD
f6 = figure('Name', 'Nichols Chart (First Harmonic) by FRD', 'NumberTitle', 'off');
nichols(sys_frd);
grid on;

toc
