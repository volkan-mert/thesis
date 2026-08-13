

w0=3.98; G = evalfr(G, 1i*w0);
%%
clear; clc; close all

R = 15;
A_i = 0.31;
w = 3.98;

varpi = R ./ (A_i * w); 

N = (4/pi) * varpi .* exp(-1i * acos((pi/2) * varpi));
