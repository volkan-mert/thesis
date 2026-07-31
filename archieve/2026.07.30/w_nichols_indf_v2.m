clear; clc; close all

%% Fig. 8 
% Negative inverse describing function technique applied to 
% the X-15 aircraft15: Nichols chart. 
% 1/N(jw; u0), negative inverse describing function and¡ 
% G(jw), linear portion of the X-15 aircraft.

%% Building Describing Function N(jw, u0) and Inverse Negative Inverse Describing Function 

% 1. Frequency ratio vector definition
n = 1000;
w = logspace(-1, 2, n);  % w = w / w_onset (0.1 to 100)
varpi = linspace(-1, 3, n);
w_onset

% 2. Preallocate vectors
A   = zeros(size(w));
phi = zeros(size(w));   % in radians

% 3. Region Masks
reg1 = (w < 1.0);
reg2 = (w >= 1.0) & (w < 1.862);
reg3 = (w >= 1.862);

% --- Region I: No saturation ---
A(reg1)   = 1.0;
phi(reg1) = 0.0;

% --- Region II: Transition (Cubic spline) ---
alpha = w(reg2);
A(reg2)   = 0.2908*alpha.^3 - 1.4396*alpha.^2 + 1.9232*alpha + 0.223;
phi(reg2) = 0.5280*alpha.^3 - 2.6213*alpha.^2 + 3.5056*alpha - 1.4171;

% --- Region III: Fully developed saturation ---
varpi = 1.0 ./ w(reg3);
A(reg3)   = (4.0 * varpi) / pi;
phi(reg3) = -acos((pi * varpi) / 2.0);


% Finding the Nichols Chart by Using FRD Object
% 
% Complex describing function N(w) = A * exp(j*phi)

N_complex       = A(reg3) .* exp(1i * phi(reg3));
w_reg3 = w(reg3);
sys_frd = frd(N_complex, w_reg3);
% sys_frd         = frd(N_complex, w);    % Describing Function
% indf_sys_frd    = -1 / sys_frd;         % Inverse Negative Describing Function
indf_sys_frd = frd(-1 ./ N_complex, w_reg3);

%% Pilot Vehicle System and their Block Diagrams

num_claw   = [5.21, -273.7855, -1425.456, -700.224];
den_claw   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw    = tf(num_claw, den_claw);

num_ldynac = [-10.524, -16.8384, -0.6209, 0];
den_ldynac = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac  = tf(num_ldynac, den_ldynac);

% G          = Gs_claw * Gs_ldynac;               % Linear OLTF q / q_c
Gs_ac          = Gs_claw*Gs_ldynac;               % Linear OLTF q / q_c

% Use freqresp to get complex response, magnitude, phase
G = freqresp(Gs_ac, w_reg3); 
G = squeeze(G);
magG = abs(G);                   % magnitude (dB)
phiG = angle(G);               % phase (rad/s)

G_complex = magG .* exp(1i * phiG);
G_sys = frd(G_complex, w_reg3);

%% Sketching the Nichols Chart
figure;
nicholsplot(indf_sys_frd, Name='Inverse Negative Describing Function');
ngrid;
title('Nichols Plot of the Negative Inverse Describing Function (by using FRD)');
legend show

hold on

% cnt = A(reg3);


nicholsplot(G_sys , Name='Aircraft Linear Transfer Function');

