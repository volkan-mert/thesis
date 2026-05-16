clear; clc; close all;
num_claw_7   = [5.21, -273.7855, -1425.456, -700.224];
den_claw_7   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw_7    = tf(num_claw_7, den_claw_7);
num_ldynac_7 = [-10.524, -16.8384, -0.6209, 0];
den_ldynac_7 = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac_7  = tf(num_ldynac_7, den_ldynac_7);
Gs_ac_7      = Gs_claw_7 * Gs_ldynac_7;
omega = logspace(-1, 2, 100);
A = ones(size(omega));
phi = zeros(size(omega));
N_frd  = frd(A .* exp(1j*phi), omega);
LN_frd = Gs_ac_7 * N_frd;

opts = nicholsoptions;
opts.PhaseMatching       = 'on';
opts.PhaseMatchingFreq   = 1;
opts.PhaseMatchingValue  = -180;
opts.PhaseWrapping       = 'on';
opts.PhaseWrappingBranch = -360;

fig = figure;
nicholsplot(LN_frd, opts); hold on;

% Manual extraction
LN_resp = squeeze(LN_frd.ResponseData);
LN_resp = LN_resp(:).';
LN_mag_dB = 20*log10(abs(LN_resp));
LN_phase_raw = rad2deg(unwrap(angle(LN_resp)));
[~, idx1] = min(abs(omega - 1));
phase_shift = -180 - LN_phase_raw(idx1);
LN_phase_matched = LN_phase_raw + phase_shift;

plot(LN_phase_matched, LN_mag_dB, 'r--', 'LineWidth', 2);
saveas(fig, 'test.png');
