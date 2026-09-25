clear; clc; close all; clear function

%% 1. Figure 8. Gain Phase Diagram for PIO Caused by a PRate-imited Servo of Ashkenas's 1964 Paper 

Kp = 1; % pilot gain
M_del_e = 1; % u_RLE: the input signal entering the rate limiter
omega_n = 2.3; % natural frequency
zeta_sp = 1.42 / omega_n / 2; % damping ratio

num = Kp*M_del_e*[1 0.82];
den = [1 1.42 omega_n^2 0];

Gs = tf(num, den) % The open-loop transfer function

%% 2. Rate limiter and input parameters

R  = 60;          % Maximum actuator rate, deg/s
% Ai = 13.68;       % Sinusoidal input amplitude, deg
Ai = Kp;

% Onset frequency
w_onset = R / Ai;     % rad/s

% Frequency vector
w = logspace(log10(0.1*w_onset), log10(10*w_onset), 2000);

% Normalized frequency
alpha = w / w_onset;

%% 3. Initialize magnitude and phase

M   = ones(size(w));
phi = zeros(size(w));

%% 4. Define the three regions

region1 = alpha <= 1;

region2 = alpha > 1 & alpha < 1.862;

region3 = alpha >= 1.862;


%% Region I: No saturation

M(region1)   = 1;
phi(region1) = 0;


%% Region II: Transition region

a = alpha(region2);

M(region2) = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.223;

phi(region2) = 0.528*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;


%% Region III: Fully developed rate saturation

varpi = w_onset ./ w(region3);

M(region3) = (4/pi) .* varpi;

phi(region3) = -acos((pi/2).*varpi);


%% 5. Complex describing function

N = M .* exp(1j*phi);

%% Quantity used for limit-cycle analysis

minus_inv_N = -1 ./ N;

%% 6. Nichols Chart

minus_inv_N_frd = frd(minus_inv_N, w);

figure
nicholsplot(Gs, minus_inv_N_frd)
grid on

title('Nichols Chart: G(j\omega) and -1/N(A,\omega)')
legend('G(j\omega)', '-1/N(A,\omega)')


%% 7. Nichols Chart Drawn Manually (SHIFTED) 
% 
% Gjw = squeeze(freqresp(Gs,w));
% 
% phase_G = rad2deg(angle(Gjw));
% gain_G = 20*log10(abs(Gjw));
% 
% phase_N = rad2deg(angle(minus_inv_N));
% gain_N = 20*log10(abs(minus_inv_N));
% 
% % Shift -1/N(A,w) to negative phase
% phase_N(phase_N > 0) = phase_N(phase_N > 0) - 360;
% 
% % Shift G(jw) to same phase range if necessary
% phase_G(phase_G > 0) = phase_G(phase_G > 0) - 360;
% 
% figure
% 
% plot(phase_G,gain_G,'LineWidth',1.5)
% hold on
% plot(phase_N,gain_N,'LineWidth',1.5)
% 
% ngrid
% 
% grid on
% 
% xlabel('Open-Loop Phase (deg)')
% ylabel('Open-Loop Gain (dB)')
% title('Nichols Chart: G(j\omega) and -1/N(A,\omega)')
% legend('G(j\omega)','-1/N(A,\omega)','Location','best')
% 
% xlim([-180 0])
%% 8. Nichols Chart (nicholsplot())

shift_dB = -17.3;

minus_inv_N_shifted = minus_inv_N * 10^(shift_dB/20);

minus_inv_N_frd = frd(reshape(minus_inv_N_shifted,1,1,[]), w);

% Nichols plot options
opt = nicholsoptions;

opt.Grid = 'on';

% Force phase into [-180, 180) deg
opt.PhaseWrapping = 'on';
opt.PhaseWrappingBranch = -180;

% Plot
figure
h = nicholsplot(Gs, minus_inv_N_frd, opt);

title('Nichols Chart: G(j\omega) and -1/N(A,\omega)')
legend('G(j\omega)', '-1/N(A,\omega)')

xlim([-180 0])