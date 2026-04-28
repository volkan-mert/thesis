clear;clc;close all
%%
tic
%

num_t_step = 1000;      % numbers of time step

X_i_max = 1;
A = X_i_max;
%
R = 0.725;        % slew rate
%
modelName = 'rle'; % Rate Limiter Element

omega = 0.1:0.1:100; % Linear increments
% omega = logspace(-1,2,100); % logarithmic increments
omega = 0.1:0.1:0.1; % Linear increments

%%

for i = 1:length(omega)

    w = omega(i);
    
    T_period = 2*pi/w;
    t_stop = T_period * 10;

    t = linspace(0, t_stop, num_t_step)'; 
    u = X_i_max*sin(w*t);        % input signal
   
    inputData = [t, u];          % [time value] format for From Workspace
    assignin('base','myInput',inputData)
    
    simOut = sim(modelName, 'StopTime', num2str(t_stop));
    
    t = simOut.tout;
    y = simOut.yout(:,1); 
    
    % idx = t > (t(end) - T_period); % Indices of the last period
    idx = t > (t_stop - T_period); % Indices of the last period
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

%% Bode Plot by FRD

figure(2);
bode(sys_frd);
hold on
xline(R/A,'r--','Label','R/A','LabelOrientation','horizontal');
grid on

%% Drawing Nyquist Plot manually

figure(3);

plot(real(H), imag(H), 'b', 'LineWidth', 1.5); 
hold on;
plot(real(H), -imag(H), 'b--', 'LineWidth', 1.5);   % mirror for -ω
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
