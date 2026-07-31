clear; clc; close all

%% --------- Linear plant : control law × aircraft dynamics ----------------
num_claw_7   = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7    = tf(num_claw_7, den_claw_7);

num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7  = tf(num_ldynac_7, den_ldynac_7);

Gs_ac_7      = Gs_claw_7 * Gs_ldynac_7;          % linear OLTF q / q_c


nopts = nicholsoptions;
nopts.PhaseMatching       = 'on';
nopts.PhaseMatchingFreq   = 1;             % rad/s
nopts.PhaseMatchingValue  = -180;          % anchor phase at -180 at w=1
nopts.PhaseWrapping       = 'on';
nopts.PhaseWrappingBranch = -360;

nicholsplot(Gs_ac_7, -1/Gs_ac_7, nopts);
legend('Aircraft Linear OLTF','DF');
grid on

ylim([-5, 15])
xlim([-180, 0])


% nicholsplot(Gs_ac_7, -1/Gs_ac_7);
% legend('Aircraft Linear OLTF','DF');