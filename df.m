%% RLE DESCRIBING FUNCTION ACCORDING TO THE ATTACHED PAPER

clear;
clc;
close all;

%% 1. Rate limiter and input parameters

R  = 60;          % Maximum actuator rate, deg/s
Ai = 13.68;       % Sinusoidal input amplitude, deg

% Onset frequency
w_onset = R / Ai;     % rad/s

% Frequency vector
w = logspace(log10(0.1*w_onset), log10(10*w_onset), 2000);

% Normalized frequency
alpha = w / w_onset;

%% 2. Initialize magnitude and phase

M   = ones(size(w));
phi = zeros(size(w));

%% 3. Define the three regions

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


%% 4. Complex describing function

N = M .* exp(1j*phi);


%% Quantity used for limit-cycle analysis

minus_inv_N = -1 ./ N;


%% 5. Important frequencies

w_crit = 1.862*w_onset;

fprintf('Input amplitude Ai       = %.4f deg\n', Ai);
fprintf('Rate limit R             = %.4f deg/s\n', R);
fprintf('Onset frequency          = %.4f rad/s\n', w_onset);
fprintf('Full-saturation boundary = %.4f rad/s\n', w_crit);


%% 6. Bode-like plot

figure(Name='Bode Plot',NumberTitle='off');

subplot(2,1,1)

semilogx(w, 20*log10(abs(N)), 'LineWidth',1.5)

grid on
hold on

xline(w_onset,'--','\omega_{onset}');
xline(w_crit,'--','1.862\omega_{onset}');

ylabel('Magnitude (dB)')
title('RLE Describing Function')


subplot(2,1,2)

semilogx(w, rad2deg(angle(N)), 'LineWidth',1.5)

grid on
hold on

xline(w_onset,'--','\omega_{onset}');
xline(w_crit,'--','1.862\omega_{onset}');

xlabel('\omega (rad/s)')
ylabel('Phase (deg)')


%% 7. Normalized frequency plot

figure(Name='Bode Plot (Normalized)',NumberTitle='off');

subplot(2,1,1)

semilogx(alpha, abs(N), 'LineWidth',1.5)

grid on
hold on

xline(1,'--');
xline(1.862,'--');

ylabel('|N(j\omega,A_i)|')
title('RLE Describing Function')


subplot(2,1,2)

semilogx(alpha,rad2deg(angle(N)), 'LineWidth',1.5)

grid on
hold on

xline(1,'--');
xline(1.862,'--');

xlabel('\omega/\omega_{onset}')
ylabel('Phase (deg)')