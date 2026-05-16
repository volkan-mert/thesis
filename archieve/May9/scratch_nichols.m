clear; clc; close all;
num_claw_7   = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7    = tf(num_claw_7, den_claw_7);
num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7  = tf(num_ldynac_7, den_ldynac_7);
Gs_ac_7      = Gs_claw_7 * Gs_ldynac_7;
omega = logspace(-1, 2, 100);

opts = nicholsoptions;
opts.PhaseMatching       = 'on';
opts.PhaseMatchingFreq   = 1;            % rad/s
opts.PhaseMatchingValue  = -180;         % anchor phase at -180° at ω = 1 rad/s
opts.PhaseWrapping       = 'on';
opts.PhaseWrappingBranch = -360;

fig = figure;
nicholsplot(Gs_ac_7, opts);
ax = gca;
lines = findobj(ax, 'Type', 'line');
disp('Lines found:');
disp(length(lines));
for i=1:length(lines)
    fprintf('Line %d points: %d\n', i, length(lines(i).XData));
end
