clear;clc;close all
%%
tic

R = 0.725;        % slew rate
%
modelName = 'rle2'; % Rate Limiter Element

%% Generating the Dominant Fourier Coefficients

% Parameters
A = 1;                     % Input Amplitude
freqs = logspace(-1, 2, 100); % 0.1 to 100 rad/s
%freqs = 0.1:0.1:100;
mag = zeros(size(freqs));
phase = zeros(size(freqs));

for i = 1:length(freqs)
    w = freqs(i);
    T_period = 2*pi/w;
    t_stop = T_period * 10; % Run for 10 periods to ensure steady state
    
    % Run Simulation
    simOut = sim(modelName, 'StopTime', num2str(t_stop));
    
    % Extract Data (assuming 'y' is the output of the nonlinearity)
    % We take only the last period to avoid transients
    t = simOut.tout;
    y = simOut.yout(:,1); 
    
    idx = t > (t_stop - T_period); % Indices of the last period
    t_last = t(idx);
    y_last = y(idx);
    
    % Perform Numerical Integration for 1st Harmonic
    % b1: Sine component, a1: Cosine component
    b1 = (2/T_period) * trapz(t_last, y_last .* sin(w * t_last));
    a1 = (2/T_period) * trapz(t_last, y_last .* cos(w * t_last));
    
    % Complex Response
    H = (b1 + 1i*a1) / A;
    
    mag(i) = abs(H); % Magnitude in dB
    phase(i) = rad2deg(angle(H)); % Phase in Degree
    G(i) = H; % Complex Response
end

%% FRD Object
sys_frd = frd(G,freqs); % frequency in rad/s

%% FIGURE 1 : LOG-MAGNITUDE & PHASE vs ω (manual)
figure('Name','Sweep Results','NumberTitle','on','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
semilogx(freqs, 20*log10(mag), 'b-o','MarkerSize',4,'LineWidth',1.5);
grid on;
xlabel('\omega  [rad/s]'); ylabel('Magnitude  [dB]');
title('Rate Limiter Describing Function – Magnitude');
xline(R/A,'r--','Label','R/A','LabelOrientation','horizontal');

ax2 = nexttile;
semilogx(freqs, phase, 'r-o','MarkerSize',4,'LineWidth',1.5);
grid on;
xlabel('\omega  [rad/s]'); ylabel('Phase  [deg]');
title('Rate Limiter Describing Function – Phase');
xline(R/A,'r--','Label','R/A','LabelOrientation','horizontal');
yline(-90,'k--','Label','-90°','LabelOrientation','horizontal');
linkaxes([ax1 ax2],'x');

%% Bode Plot by FRD

figure(2);
bode(sys_frd);
hold on
xline(R/A,'r--','Label','R/A','LabelOrientation','horizontal');
grid on

%% Drawing Nyquist Plot manually

figure(3);

plot(real(G), imag(G), 'b', 'LineWidth', 1.5); 
hold on;
plot(real(G), -imag(G), 'b--', 'LineWidth', 1.5);   % mirror for -ω
hold on
plot(-1, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 2); % critical point
grid on; 
axis equal;
xlabel('Real Axis'); ylabel('Imaginary Axis');
title('Nyquist Diagram');

%% Drawing Nyquist Plot by using FRD

figure(4);

nyquist(sys_frd); % Nyquist Plot
grid on;

%% Drawing Nichols Chart manually

figure(5);

plot(phase, mag, 'b', 'LineWidth', 1.5);
hold on
plot(-180, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);  % critical point
xlabel('Open-Loop Phase (deg)');
ylabel('Open-Loop Gain (dB)');
title('Nichols Chart');
grid on;
xlim('auto');
ylim('auto');
% xlim([-270 -90]);
% ylim([-40 20]);

%% Drawing Nichols Chart by using FRD 

figure(6);

nichols(sys_frd);
grid on

%
toc