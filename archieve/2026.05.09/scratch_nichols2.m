clear; clc; close all;
num = [1]; den = [1 1];
G = tf(num, den);
omega = logspace(1, 2, 10);
G_frd = frd(G, omega);

opts = nicholsoptions;
opts.PhaseMatching       = 'on';
opts.PhaseMatchingFreq   = 1;
opts.PhaseMatchingValue  = -180;

try
    nicholsplot(G_frd, opts);
    disp('Success!');
catch ME
    disp(['Error: ', ME.message]);
end
