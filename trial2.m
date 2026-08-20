
clear;
clc;
close all;

%%

numGc = [5.21, -273.7855, -1425.2, -700.1952];
denGc  = [1, 21.3594, 545.5538, 605.6621, 0];

Gac = tf(numGac,denGac);

%% 2. Control Law

numGac = [-10.5240,-16.8384,-0.6247,0];
denGac = [1, 2.3473, -5.3061, -0.1836, -0.0418];

Gc = tf(numGc,denGc);


num_ctrl = 5.21 * [1, -52.55, -273.6, -134.4];
den_ctrl = [1, 21.36, 545.6, 605.7, 0];

sys_ctrl = tf(num_ctrl, den_ctrl);

%% 3. Longitudinal Dynamics of the Aircraft
num_dyn = -10.524 * [1, 1.6, 0.059, 0];
den_dyn = [1, 2.35, -5.31, 0.184, -0.041];

sys_dyn = tf(num_dyn, den_dyn);

%%

isequal(Gac, sys_dyn)
isequal(Gc,  sys_ctrl)


%% 
 
all(abs(Gac.Numerator{1} - sys_dyn.Numerator{1}) < 1e-10) && all(abs(Gac.Denominator{1} - sys_dyn.Denominator{1}) < 1e-10)

%%

isequal(sort(pole(Gac)), sort(pole(sys_dyn)))
isequal(sort(zero(Gc)), sort(zero(sys_ctrl)))
