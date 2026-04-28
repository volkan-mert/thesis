clear; clc; close all
%%

modelName = 'rle2';         % simulink model of sine wave
A = 1;                      % amplitude of sine wave
omega = logspace(-1,2,100); % angular frequency
mag = zeros(size(omega));   % placeholder of magnitude
phase = zeros(size(omega)); % placeholder of phase
R = 0.725;                  % rising slew rate

for i = 1:length(omega)

    w = omega(i); % angular frequency in a cell
    
    T = 2*pi/w; % period in a cell
    
    t_stop = T*10; % run for 10 periods to ensure the steady state

    out = sim(modelName, 'StopTime',num2str(t_stop)); % runnig the simulation

    t = out.tout; % time series
    y = out.yout; % result of the sine wave through rate limiter

    idx = t > (t(end)- T); % if the value of t is greater than or equal to the period takes it as the present time

    t_last = t(idx);
    y_last = y(idx);

    % Perform Numerical Integration for 1st Harmonic
    
    b1 = (2/T)*trapz(t_last, y_last.*sin(w*t_last)); % sine component
    a1 = (2/T)*trapz(t_last, y_last.*cos(w*t_last)); % cosine component

    H = (b1 + 1i*a1) / A; % complex response

    mag(i) = abs(H); % magnitude in dB
    phase(i) = rad2deg(angle(H)); % phase in degree
    G(i) = H; % complex response
    
end

%% Frequency Response Data-Model

sys_frd = frd(G, omega); % frequency response data-model

%% Drawing Bode Diagram manually

figure(1);

subplot(2,1,1);
semilogx(omega, 20*log10(mag), 'b-o', 'LineWidth', 1.5);
grid on; 
ylabel('Magnitude (dB)'); 
title('Nonlinear Bode Plot (First Harmonic)');

subplot(2,1,2);
semilogx(omega, phase, 'r-o', 'LineWidth', 1.5);
grid on; 
ylabel('Phase (deg)'); 
xlabel('Frequency (rad/s)');

%% Drawing Bode Diagram by using FRD

figure(2);

bode(sys_frd)
xline(R/A,'DisplayName','w_onset');
grid on;


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


