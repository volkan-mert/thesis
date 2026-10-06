% Describing Function Analysis - Limit Cycle via Nested Golden-Section Search
clear; close all; clc

% Parameters
mu  = 1;
tol = 1e-8;                          % golden-section tolerance

% The two curves
G   = @(w)    mu ./ ((1 - w.^2) - 1j*mu*w);
niN = @(A, w) 1j * 4 ./ (A.^2 .* w);

% Distance between the two curves at the same frequency w
dist = @(w, A) abs(G(w) - niN(A, w));

% Search brackets (read roughly from the Nyquist plot)
w_lo = 0.5;  w_hi = 1.4;
A_lo = 1;    A_hi = 5;

% Inner search: best A for a given w, returns the smallest distance
bestA = @(w) golden(@(A) dist(w, A), A_lo, A_hi, tol);
phi   = @(w) dist(w, bestA(w));

% Outer search: frequency where the smallest distance is minimum
w_star = golden(phi, w_lo, w_hi, tol);
A_star = bestA(w_star);
Gs     = G(w_star);                  % intersection point

fprintf('Limit cycle frequency  w* = %.8f rad/s\n', w_star)
fprintf('Limit cycle amplitude  A* = %.8f\n', A_star)
fprintf('Intersection point        = %.8f + %.8fj\n', real(Gs), imag(Gs))
fprintf('Distance at solution      = %.2e\n', dist(w_star, A_star))

% Figure 1: full Nyquist plot
% Frequency grid that contains w* exactly
w = unique([0.01:0.1:5, w_star]);

sys_G   = frd(G(w), w);
sys_niN = frd(niN(A_star, w), w);

figure
nyquist(sys_G, 'b')
hold on
grid on
nyquist(sys_niN, 'r*')
plot(real(Gs), imag(Gs), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5)
hold off
ylim([-10 10])
title('The Occurrence of the Limit Cycle')
legend('G(jw)', sprintf('-1/N, A = %.4f', A_star), 'Intersection')

% Figure 2: zoomed view around the intersection
% Fine frequency grid near w* so the curves are accurate inside the small window
w_zoom = unique([w_star + (-0.002:0.0001:0.002), w_star]);

sys_G_zoom   = frd(G(w_zoom), w_zoom);
sys_niN_zoom = frd(niN(A_star, w_zoom), w_zoom);

figure
nyquist(sys_G_zoom, 'b')
hold on
grid on
nyquist(sys_niN_zoom, 'r*')
plot(real(Gs), imag(Gs), 'ko', 'MarkerSize', 12, 'LineWidth', 1.5)
hold off
xlim([real(Gs) - 0.001, real(Gs) + 0.001])
ylim([imag(Gs) - 0.001, imag(Gs) + 0.001])
title(sprintf('Zoom at Intersection: w* = %.4f rad/s, A* = %.4f', w_star, A_star))
legend('G(jw)', sprintf('-1/N, A = %.4f', A_star), 'Intersection')

% Golden-section search (local function)
function xmin = golden(f, a, b, tol)
r  = (sqrt(5) - 1) / 2;          % about 0.618
x1 = b - r*(b - a);
x2 = a + r*(b - a);
f1 = f(x1);
f2 = f(x2);
while (b - a) > tol
    if f1 < f2                   % minimum is in [a, x2]
        b  = x2;
        x2 = x1;  f2 = f1;
        x1 = b - r*(b - a);  f1 = f(x1);
    else                         % minimum is in [x1, b]
        a  = x1;
        x1 = x2;  f1 = f2;
        x2 = a + r*(b - a);  f2 = f(x2);
    end
end
xmin = (a + b) / 2;
end