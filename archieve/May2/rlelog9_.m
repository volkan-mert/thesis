clear;clc;close all
%%
tic
%
t_stop = 10*2*pi;      % simulation stop time
num_t_step = 1e4;      % numbers of time step
t = linspace(0, t_stop, num_t_step)';      % time
X_i_max = 1;
A = X_i_max;
%
R = 0.725;        % slew rate
%
modelName1 = 'rle9sinewave'; % sine wave generation
modelName2 = 'rle9'; % Rate Limiter Element

load_system(modelName1); % Loads the model into memory without opening the window
load_system(modelName2); % Loads the model into memory without opening the window

% Show Simulink Model and the Scope

% open_system(modelName);

    % scopes = find_system(modelName2, 'BlockType', 'Scope');
    % for sc = 1:length(scopes)
    %     open_system(scopes{sc});
    % end

% %% Blocks of Simulink Model
% 
% % Aircraft Block of the aircraft transfer function of X-15 PIO during flare
% % (2000, F. Amato*, R. Iervolino*, M. Pandit†, S.Scala‡, L.Verde, Analysis of Pilot-in-the-Loop Oscillations Due to Position and Rate Saturations, CDC Paper)
% 
% num_ac = [3.476 3.1708072 0.0896237936]; % numerator of the transfer function of the Aircraft
% 
% den_ac = [1 0.8608 5.3159942 0.108928 0.0529]; % denominator of the transfer function of the Aircraft
% 
% Gs_ac = tf(num_ac, den_ac); % transfer function of the aircraft
% 
% % Pilot Gain Block
% Kp = 3.5; % pilot gain (corresponding to 20 degrees of phase margin of the linearised model)
Kp = 1; % pilot gain

%%

cnt = 0;

% omega = 0.1:0.1:100; % Linear increments
omega = logspace(-1,2,100); % logarithmic increments
% omega = 0.1:0.1:0.1; % Linear increments
u = zeros(num_t_step,length(omega));

for i = 1:length(omega)
    w = omega(i);
    % u(:, i) = X_i_max*sin(w*t);        % input signal
    out1 = sim(modelName1, 'StopTime', num2str(t_stop));
    % u(:, i) = out1.myinput.Data;
    u(:, i) = interp1(out1.myinput.Time, out1.myinput.Data, t, 'linear', 'extrap');
end

% inputData = [t, u];          % [time value] format for From Workspace
% assignin('base','myInput',inputData);

inputDataTS = timeseries(u, t);     % u is num_t_step × length(omega), t is num_t_step × 1
assignin('base','myInput',inputDataTS);

if any(isnan(u(:))) || any(isinf(u(:)))   % must be 0
    disp('Information: Input Signal Data contain Inf or NaN values!!! Change or fix the variable''inputDataTS''.');
else
    disp('Information: Input Signal Data don''t contain Inf or NaN values.');
end

% In the model, set From Workspace block Data = myInput
out2 = sim(modelName2);
% After sim finishes, read output variable from To Workspace, e.g.:
outData = out2.yout;           % if To Workspace used default name
t = out2.tout;

%%

for i = 1:length(omega)

    w = omega(i);
    y = outData(:, i);
    T_period = 2*pi/w;
    t_stop = T_period * 10;
    
    % idx = t > (t(end) - T_period); % Indices of the last period
    idx = t > (t(end) - T_period); % Indices of the last period
    t_last = t(idx);
    y_last = y(idx);

    % Perform Numerical Integration for 1st Harmonic
    % b1: Sine component, a1: Cosine component
    b1 = (2/T_period) * trapz(t_last, y_last .* sin(w * t_last));
    a1 = (2/T_period) * trapz(t_last, y_last .* cos(w * t_last));

    % Complex Response
    H(i) = (b1 + 1i*a1) / A;

    mag(i) = abs(H(i));
    phase(i) = rad2deg(angle(H(i)));

end

%% Frequency Response Data-model

sys_frd = frd(H, omega);

%% Bode Plot
figure(1);
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

%%

pick_cursor_magnitude_olop = 2.0629; % omega_OLOP Magnitude in dB
pick_cursor_phase_olop = -136.6216; % omega_OLOP Phase in degrees

%%

% --- Clean the describing-function data ---
% (1) Force exact unity below the activation onset (no rate limit there)
inactive          = (omega * A) <= R;
H(inactive)       = abs(H(inactive));     % zero out the numerical imag noise

% (2) Unwrap the angle so it's monotonic across omega
phi_unw           = unwrap(angle(H));
H_clean           = abs(H) .* exp(1i*phi_unw);

sys_frd           = frd(H_clean, omega);

% --- Plot with phase matching anchored where BOTH curves are well-defined ---
figure(2);
h1 = nicholsplot(sys_frd);   p1 = getoptions(h1);
p1.PhaseMatching       = 'on';
p1.PhaseMatchingFreq   = 0.1;       % low freq: DF phase ≈ 0
p1.PhaseMatchingValue  = 0;
p1.PhaseWrapping       = 'on';
p1.PhaseWrappingBranch = -360;
setoptions(h1,p1);

hold on

h2 = nicholsplot(Gs_ac_7);   p2 = getoptions(h2);
p2.PhaseMatching       = 'on';
p2.PhaseMatchingFreq   = 0.1;       % SAME anchor frequency
p2.PhaseMatchingValue  = rad2deg(angle(evalfr(Gs_ac_7, 1i*0.1)));   % its true phase there
p2.PhaseWrapping       = 'on';
p2.PhaseWrappingBranch = -360;
setoptions(h2,p2);


%%
plot(pick_cursor_phase_olop, pick_cursor_magnitude_olop,'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
hold on
text(pick_cursor_phase_olop + 5, pick_cursor_magnitude_olop, '\omega_{OLOP} = 5.1', 'Color', 'red', 'FontWeight', 'bold', 'FontSize', 12);
legend('The 1st Harmonics','CLAW and Longitudinal Dynamics of the Aircraft');
grid on


