
clear; clc; close all;

%% 1. Forward Gain
Kf = 13.68;

%% 2. Control Law
num_ctrl = 5.21 * [1, -52.55, -273.6, -134.4];
den_ctrl = [1, 21.36, 545.6, 605.7, 0];

sys_ctrl = tf(num_ctrl, den_ctrl);

%% 3. Longitudinal Dynamics of the Aircraft
num_dyn = -10.524 * [1, 1.6, 0.059, 0];
den_dyn = [1, 2.35, -5.31, 0.184, -0.041];

sys_dyn = tf(num_dyn, den_dyn);

%% 4. Complete Linear System
sys = Kf * sys_ctrl * sys_dyn;


nyquist(sys)