%% Nichols Chart: Negative Inverse Describing Function Analysis
% Replicating Fig. 8 from Duda (1997) for custom Simulink TF

clear; clc; close all;

%% 1. Define the Linear Transfer Function G(s)
% Convolving CLAW and Aircraft transfer functions as provided
num = conv([5.21, -273.7855, -1425.2403, -700.1952], [-10.524, -16.8384, -0.6247]);
den = conv([1, 21.3594, 545.5538, 605.6621], [1, 2.3473, -5.3061, -0.1836, -0.0418]);
G = tf(num, den);

nicholsplot(G, -1/G)
legend('Linear A/C TF','Describing Function')
grid on