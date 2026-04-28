clear all, clc; close all;

omega_onset = 1;
omega = (pi/2)*omega_onset:0.01:10;
N=(4/pi)*(omega_onset./omega).*exp(-1j*acos((pi/2)*(omega_onset./omega)));
Ng=20*log(abs(N));
plot(omega,Ng)
