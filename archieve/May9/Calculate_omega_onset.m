clear;clc;close all


%% Example : Draw Fig. 16 Nichols chart: open-loop frequency response q∕qc 
% 2025, Rodrigues et al, CEAS Aeronautical Journal (2026), 
% "Conditional integrator sliding mode control to reduce susceptibility 
% to pilot‑induced oscillations"


num_claw_7 = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7 = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7 = tf(num_claw_7, den_claw_7);


num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7 = tf(num_ldynac_7, den_ldynac_7);


F_uc_2_urle = feedback(Gs_claw_7,Gs_ldynac_7);


% Double check the value of the closed-loop onset frequency
R = deg2rad(60);
w=5.077;
Ain = 13.68/180*pi;  % since pilot input is 100% and kf=13.68
rate = Ain*abs(evalfr(F_uc_2_urle,1j*w))*w;


% Find the closed-loop onset frequency
f = @(w) (R - Ain*abs(evalfr(F_uc_2_urle,1j*w))*w)^2;
A = []; b = []; Aeq = []; beq = [];
lb = 0; ub = 100;
w0 = 0;
[w_opt, fval] = fmincon(f, w0, A, b, Aeq, beq, lb, ub)
