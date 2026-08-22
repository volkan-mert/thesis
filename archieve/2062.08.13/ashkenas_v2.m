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

shift_dB = -18;

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

%% 9. Figure 5.24 - G(jw) and -1/N(A,w)
% For this figure, frequency is fixed for each -1/N curve.
% Amplitude A changes along each curve.

w_G = logspace(-2,2,2000);
[Re_G, Im_G] = nyquist(Gs,w_G);

Re_G = squeeze(Re_G);
Im_G = squeeze(Im_G);

w_fixed = [1 2 3 4];          % Fixed frequencies, rad/s
A_values = linspace(0.01,90,1000);

figure
plot(Re_G,Im_G,'k','LineWidth',1.5)
hold on

for k = 1:length(w_fixed)

    w0 = w_fixed(k);

    alpha_24 = A_values*w0/R;

    M_24 = ones(size(alpha_24));
    phi_24 = zeros(size(alpha_24));

    region1 = alpha_24 <= 1;
    region2 = alpha_24 > 1 & alpha_24 < 1.862;
    region3 = alpha_24 >= 1.862;

    % Region I
    M_24(region1) = 1;
    phi_24(region1) = 0;

    % Region II
    a = alpha_24(region2);

    M_24(region2) = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.223;

    phi_24(region2) = 0.528*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;

    % Region III
    varpi = 1./alpha_24(region3);

    M_24(region3) = (4/pi).*varpi;

    phi_24(region3) = -acos((pi/2).*varpi);

    % Describing function
    N_24 = M_24.*exp(1j*phi_24);

    % -1/N(A,w)
    minus_inv_N = -1./N_24;

    plot(real(minus_inv_N),imag(minus_inv_N),'LineWidth',1.2)

end

% Critical point
plot(-1,0,'kx','MarkerSize',8,'LineWidth',1.5)

text(-1,0,'  -1')

grid on
axis equal

xlabel('Real Axis')
ylabel('Imaginary Axis')

title('Figure 5.24: G(j\omega) and -1/N(A,\omega)')

legend('G(j\omega)', '\omega = 1', '\omega = 2', '\omega = 3', '\omega = 4', '-1', 'Location','best')


%% 10. Figure 5.25 - G(jw)N(A,w)
% For this figure, A is fixed for each curve.
% Frequency changes from a small value to a large value.

A_fixed = [5 10 20 40];       % Fixed amplitudes, deg

w_25 = logspace(-2,3,3000);    % Frequency vector, rad/s

% Calculate G(jw)
Gjw = squeeze(freqresp(Gs,w_25));

% Convert G(jw) to a row vector
Gjw = reshape(Gjw,1,[]);

L = cell(1,length(A_fixed));

for k = 1:length(A_fixed)

    A0 = A_fixed(k);

    alpha_25 = A0*w_25/R;

    M_25 = ones(size(alpha_25));
    phi_25 = zeros(size(alpha_25));

    region1 = alpha_25 <= 1;
    region2 = alpha_25 > 1 & alpha_25 < 1.862;
    region3 = alpha_25 >= 1.862;

    % Region I
    M_25(region1) = 1;
    phi_25(region1) = 0;

    % Region II
    a = alpha_25(region2);

    M_25(region2) = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.223;

    phi_25(region2) = 0.528*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;

    % Region III
    varpi = 1./alpha_25(region3);

    M_25(region3) = (4/pi).*varpi;

    phi_25(region3) = -acos((pi/2).*varpi);

    % Describing function
    N_25 = M_25.*exp(1j*phi_25);

    % Convert N(A,w) to a row vector
    N_25 = reshape(N_25,1,[]);

    % Calculate G(jw)N(A,w)
    GN = Gjw.*N_25;

    % Convert the result to FRD form for nyquistplot()
    L{k} = frd(reshape(GN,1,1,[]),w_25);

end

% Nyquist plot options
opt = nyquistoptions;

opt.Grid = 'on';

opt.ShowFullContour = 'off';

figure

nyquistplot(L{:},opt)

hold on

% Critical point
plot(-1,0,'kx','MarkerSize',8,'LineWidth',1.5)

text(-1,0,'  -1')

title('Figure 5.25: G(j\omega)N(A,\omega)')

legend('A = 5 deg', 'A = 10 deg', 'A = 20 deg', 'A = 40 deg', 'Location','best')