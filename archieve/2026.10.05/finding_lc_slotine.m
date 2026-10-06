% Van der Pol Oscillator - Limit Cycle Analysis
% Compares the numerical limit cycle with the describing function (DF)
% prediction (Slotine & Li, Applied Nonlinear Control, Fig. 5.24-5.25)
% Requires the Control System Toolbox (tf, frd, nyquistplot)

clear; clc; close all;

% Use LaTeX for legends and tick labels.
% Titles, axis labels and text get LaTeX through ltx{:} below.
% (A global defaultTextInterpreter = 'latex' breaks the Nyquist grid
%  labels such as '-90 deg', so it is NOT set here.)
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
ltx = {'Interpreter', 'latex'};

%% Parameters
mu = 1;            % Van der Pol parameter
tspan = [0 60];    % simulation time (s)
x0 = [0.5; 0];     % initial conditions [x; dx/dt]

%% Simulation
% Equation: x'' + mu*(x^2 - 1)*x' + x = 0
% States:   x1 = x,  x2 = dx/dt
f = @(t,x) [x(2); mu*(1 - x(1)^2)*x(2) - x(1)];
[t, x] = ode45(f, tspan, x0);

%% Limit cycle amplitude and frequency from simulation
% Use only the data after t = 40 s, when the transient has died out
idx = t >= 40;
t2 = t(idx);
x2 = x(idx,1);
v2 = x(idx,2);

% Amplitude = half of the peak-to-peak value
A_num = (max(x2) - min(x2))/2;

% Period = time between two upward zero crossings
zc = find(x2(1:end-1) < 0 & x2(2:end) >= 0);
T_num = mean(diff(t2(zc)));
w_num = 2*pi/T_num;

fprintf('Simulation: A = %.3f, w = %.3f rad/s, T = %.3f s\n', A_num, w_num, T_num);
fprintf('DF method:  A = 2.000, w = 1.000 rad/s\n');

%% Time response
figure
plot(t, x(:,1), 'LineWidth', 1.5)
grid on
xlabel('$t$ (s)', ltx{:}), ylabel('$x(t)$', ltx{:})
title('Van der Pol Oscillator -- Time Response', ltx{:})

%% Phase portrait
th = linspace(0, 2*pi, 200);

figure
plot(x(:,1), x(:,2)), hold on
plot(x2, v2, 'r', 'LineWidth', 2)
plot(2*cos(th), 2*sin(th), 'k--', 'LineWidth', 1.5)   % DF circle, A = 2
grid on, axis equal
xlabel('$x$', ltx{:}), ylabel('$\dot{x}$', ltx{:})
title('Van der Pol Oscillator -- Phase Portrait', ltx{:})
legend('Trajectory', 'Limit cycle (simulation)', 'Limit cycle (DF, $A = 2$)')

%% Describing function analysis
% Linear part:      G(s) = 1/(s^2 - mu*s + 1)
% Describing func.: N(A,w) = j*mu*w*A^2/4   ->   N(A,s) = (mu*A^2/4)*s
% Limit cycle when: G(jw)*N(A,w) = -1
s = tf('s');
G = 1/(s^2 - mu*s + 1);
w = linspace(0.05, 3, 1000);      % frequency range (rad/s)

% Common Nyquist plot options (LaTeX labels, w > 0 branch only)
opts = nyquistoptions;
opts.ShowFullContour = 'off';
opts.Grid = 'on';
opts.XLabel.String = '$\mathrm{Re}$';
opts.YLabel.String = '$\mathrm{Im}$';
opts.XLabel.Interpreter = 'latex';
opts.YLabel.Interpreter = 'latex';
opts.Title.Interpreter = 'latex';

%% Fig. 5.24: G(jw) and -1/N(A,w) for fixed frequencies
% For a fixed w, -1/N depends only on A, so it is not a transfer function.
% It is stored as an frd object whose "frequency" vector holds the A values.
A = linspace(0.5, 4, 400);
wlist = [0.5 0.8 1 1.5 2];

Gf = frd(G, w);
M = cell(1, length(wlist));
for k = 1:length(wlist)
    N = 1i*mu*wlist(k)*A.^2/4;
    M{k} = frd(-1./N, A);
end

opts.Title.String = '$G(j\omega)$ and $-1/N(A,\omega)$ for fixed $\omega$ (Slotine''s Book Fig. 5.24)';
opts.XLim = {[-1 1.5]};       % zoom in on the intersection region
opts.YLim = {[-0.5 4]};

figure
nyquistplot(Gf, '-', M{1}, '--', M{2}, '--', M{3}, '--', M{4}, '--', M{5}, '--', opts)
hold on
plot(0, 1/mu, 'ko', 'MarkerFaceColor', 'k')   % solution A = 2, w = 1
legend('$G(j\omega)$', '$\omega = 0.5$', '$\omega = 0.8$', '$\omega = 1$', ...
       '$\omega = 1.5$', '$\omega = 2$', '$A = 2,\ \omega = 1$')

% Annotations
% Mark G(jw) at the same frequencies used for the -1/N curves
Gw = squeeze(freqresp(G, wlist));
plot(real(Gw), imag(Gw), 'ks', 'HandleVisibility', 'off')
for k = [1 2 4 5]             % w = 1 is labelled separately below
    text(real(Gw(k)) + 0.05, imag(Gw(k)), sprintf('$\\omega = %.1f$', wlist(k)), ltx{:})
end
text(0.08, 1/mu, '$\leftarrow$ Limit cycle: $A = 2$, $\omega = 1$', ltx{:})
text(0.08, 3, '$\downarrow$ $A$ increases along $-1/N$', ltx{:})
text(0.95, 0.3, '$\nwarrow$ $\omega$ increases along $G(j\omega)$', ltx{:})
text(-0.95, -0.3, 'Only at $\omega = 1$ does $G(j\omega)$ reach the $-1/N$ line', ltx{:})

%% Fig. 5.25: G(jw)*N(A,w) for fixed amplitudes
% For a fixed A, G(s)*N(A,s) is a normal transfer function
Alist = [1 1.5 2 2.5 3];

L = cell(1, length(Alist));
for k = 1:length(Alist)
    L{k} = G * (mu*Alist(k)^2/4) * s;
end

opts.Title.String = '$G(j\omega)N(A,\omega)$ for fixed $A$ (Slotine''s Book Fig. 5.25)';
opts.XLim = {[-2.6 0.4]};
opts.YLim = {[-1.4 1.4]};

figure
nyquistplot(L{:}, w, opts)
hold on
plot(-1, 0, 'ko', 'MarkerFaceColor', 'k')     % critical point
legend('$A = 1$', '$A = 1.5$', '$A = 2$', '$A = 2.5$', '$A = 3$', '$(-1,\,0)$')

% Annotations
% At w = 1:  G(j1) = j/mu  and  N = j*mu*A^2/4,  so  G*N = -A^2/4 (real axis)
plot(-Alist.^2/4, zeros(size(Alist)), 'k.', 'MarkerSize', 14, 'HandleVisibility', 'off')
text(-1, -0.2, '$\uparrow$ $A = 2$ passes through $(-1,0)$: limit cycle', ...
     'HorizontalAlignment', 'center', ltx{:})
text(-1.125, 1.25, '$\leftarrow$ $\omega$ increases', 'HorizontalAlignment', 'center', ltx{:})
text(-2.55, 1.25, '$\bullet$ $\omega = 1$: $GN = -A^2/4$', ltx{:})
text(-2.55, -1.25, 'Larger $A$ $\Rightarrow$ larger circle', ltx{:})

%% Fig. 5.26: Stability of the limit cycle
% From Fig. 5.25 the DF limit cycle is A = 2, w = 1 rad/s.
% With N(A,w), the system behaves like the linear equation
%     x'' + c_eff*x' + x = 0,    c_eff = mu*(A^2/4 - 1)
% c_eff < 0 : negative damping -> amplitude grows
% c_eff > 0 : positive damping -> amplitude decays
% Stable limit cycle: A grows when below it and decays when above it.
A_lc = 2;  w_lc = 1;
dA = 0.1;
A_low  = A_lc - dA;
A_high = A_lc + dA;

c_low  = mu*(A_low^2/4  - 1);
c_high = mu*(A_high^2/4 - 1);

if c_low < 0 && c_high > 0
    stability = 'stable';
else
    stability = 'unstable';
end

fprintf('\nStability check, c_eff = mu*(A^2/4 - 1):\n');
fprintf('A = %.2f (below): c_eff = %+.3f\n', A_low, c_low);
fprintf('A = %.2f (above): c_eff = %+.3f\n', A_high, c_high);
fprintf('The limit cycle is %s.\n', stability);

% -1/N(A, w = 1) for A around the limit cycle (stored as frd, as in Fig. 5.24)
A_s = linspace(A_lc - 0.5, A_lc + 0.5, 300);
Ms  = frd(-1./(1i*mu*w_lc*A_s.^2/4), A_s);

% -1/N at the three amplitudes
Mp = -1./(1i*mu*w_lc*[A_low A_lc A_high].^2/4);

opts.Title.String = ['Limit cycle stability: ' stability ' (Slotine''s Book Fig. 5.26)'];
opts.XLim = {[-0.8 1.2]};
opts.YLim = {[0.4 1.6]};

figure
nyquistplot(Gf, '-', Ms, '--', opts)
hold on
plot(real(Mp(1)), imag(Mp(1)), 'bo', 'MarkerFaceColor', 'b')
plot(real(Mp(2)), imag(Mp(2)), 'ko', 'MarkerFaceColor', 'k')
plot(real(Mp(3)), imag(Mp(3)), 'ro', 'MarkerFaceColor', 'r')
legend('$G(j\omega)$', '$-1/N(A,\ \omega = 1)$', ...
       sprintf('$A = %.1f$ (below)', A_low), '$A = 2$ (limit cycle)', ...
       sprintf('$A = %.1f$ (above)', A_high))

% Annotations: both neighbours move back toward the limit cycle point
text(0.05, imag(Mp(1)), '$\downarrow$ $c_{\mathrm{eff}} < 0$: $A$ grows', ltx{:})
text(0.05, imag(Mp(3)), '$\uparrow$ $c_{\mathrm{eff}} > 0$: $A$ decays', ltx{:})
text(-0.05, imag(Mp(2)), 'Limit cycle $A = 2$, $\omega = 1$ $\rightarrow$', 'HorizontalAlignment', 'right', ltx{:})