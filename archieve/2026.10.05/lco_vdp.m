% Describing Function Analysis - Limit Cycle Prediction (Nyquist Plot)
clear; close all; clc
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')

% Parameters
mu = 1;                 % system parameter
w  = (1:500)/100;       % frequency vector [rad/s], 0.01 ... 5, contains exactly w = 1
A  = (1:500)/100;       % amplitude values,        0.01 ... 5, contains exactly A = 2

% Linear part G(jw)
G     = mu ./ ((1 - w.^2) - 1j*mu*w);
sys_G = frd(G, w);

figure
nyquist(sys_G, 'b')
hold on
grid on

% Negative inverse describing function -1/N(A,w) for each amplitude
allZero = true;         % assume all real parts are zero until proven otherwise

for k = 1:length(A)
    niN     = 1j * 4 ./ (A(k)^2 * w);
    sys_niN = frd(niN, w);
    nyquist(sys_niN, 'r*')

    % Check the real part of -1/N
    if any(real(niN) ~= 0)
        allZero = false;
    end
end

hold off
title('The Occurrence of the Limit Cycle')
legend('$G(j\omega)$', '$-1/N(A_n, \omega)$', 'Location', 'bestoutside', 'AutoUpdate', 'off')

% Zoom in around the expected intersection point (0, j1)
xlim([-0.005, 0.005])
ylim([0.995, 1.005])

% Print a note in the Command Window
if allZero
    disp('NOTE:')
    disp('The Negative Inverse Describing Function curves do not have any real parts,')
    disp('and Re{-1/N(A, w)} = 0 for all values, which means that all of these curves overlap.')
    disp(' ')
end

% Numerical check of the intersection at w = 1 rad/s and A = 2
iw = find(w == 1);      % index of w = 1 rad/s
iA = find(A == 2);      % index of A = 2

G_1   = G(iw);
niN_1 = 1j * 4 / (A(iA)^2 * w(iw));

disp('Intersection check at w = 1 rad/s:')
fprintf('  G(j1)           = %.4f + j%.4f\n', real(G_1),   imag(G_1))
fprintf('  -1/N(A = 2, 1)  = %.4f + j%.4f\n', real(niN_1), imag(niN_1))
fprintf('  Difference      = %.2e\n', abs(G_1 - niN_1))
disp(' ')
fprintf('Predicted limit cycle: A = %.2f, w = %.2f rad/s\n', A(iA), w(iw))