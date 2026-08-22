clear; clc; close all; clear function

%% 1. Figure 8. Gain Phase Diagram for PIO Caused by a PRate-imited Servo of Ashkenas's 1964 Paper 

Kp = 1; % pilot gain
M_del_e = 1; % u_RLE: the input signal entering the rate limiter
omega_n = 2.3; % natural frequency
zeta_sp = 1.42 / omega_n / 2; % damping ratio

num = Kp*M_del_e*[1 0.82];
den = [1 1.42 omega_n^2 0];

Gs = tf(num, den) % The open-loop transfer function