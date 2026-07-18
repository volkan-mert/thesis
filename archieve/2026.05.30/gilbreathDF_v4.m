%% Reference
%  Gregory P. Gilbreath, Captain, USAF, Jan 26, 2001
%
%  "Prediction of Pilot-Induced Oscillations (PIO)
%   Due to Actuator Rate Limiting Using the
%   Open-Loop Onset Point (OLOP) Criterion"
%
%  Thesis  AFIT/GAE/ENY/01M-02
%
%  Closed-loop describing function for the aircraft model from
%  AIAA-95-3204-CP (Duda / Hanke / Gilbreath, OLOP criterion for the
%  rate-limited PIO problem).
%
%  Fixes vs. original:
%    1. eqs() referenced dfuncmag / dfuncphase, which were never defined.
%       It now calls dfunction() once and uses both outputs (keeps
%       magnitude and phase on the same spline branch).
%    2. bode() outputs at a single frequency are squeeze()d to scalars.
%    3. Dead "phi2(1);" statement removed.
%    4. Solver discontinuity is reported once (not overwritten each fail).
%    5. Arrays pre-allocated.
%    6. Nichols plot labelled.
%    7. Figure(1) is plotted with nicholsplot instead of plot.
%    8. The graphs of figure(1) have been aligned to 0 dB / -180 deg.
%    9. Linear OLOP onset frequency: onset_freq is now defined by
%       w * |P_RLE(jw)| = R   (the qco-independent definition used in
%       the OLOP framework). The fminsearch-failure point is kept as
%       solver_jump_freq for solver diagnostics.

clear; clc; close all
global w R K magGc magGac phGc phGac dtr qco z

n = 10000;                                    % frequency-grid resolution

% System input frequency range, command amplitude, rate limit
w   = logspace(-1, 2, n);
qco = 1.1;
R   = 60;

% System numbers
K   = 13.68;
Gc  = tf( 5.21 * conv([1 -57.36], conv([1 4.26], [1 0.55])), ...
          conv([1 2*0.442*22.85 22.85^2], conv([1 0], [1 1.16])) );
Gac = tf(-10.524 * conv([1 1.562], conv([1 0.038], [1 0])), ...
          conv([1 2*0.212*0.088 0.088^2], conv([1 3.75], [1 -1.44])) );

% Linear closed-loop response from qco to the rate-limiter input:
%   pcl = K*Gc / (1 + Gc*Gac)  ==  P_RLE in the OLOP nomenclature
pcl = K*Gc / (1 + Gc*Gac);
dtr = pi/180;

% Pre-allocate
N      = numel(w);                  % number of frequency points
deltao = zeros(1, N+1);             % surface-command amplitude at RLE input
phi2   = zeros(1, N+1);             % phase of the surface command [deg]
fval   = zeros(1, N);               % harmonic-balance residual
NoMag  = zeros(1, N);               % open-loop DF magnitude [dB]
NoPh   = zeros(1, N);               % open-loop DF phase [rad]

% Linear initial guess (low-frequency point)
[magpcl, p0] = bode(pcl, w(1));
deltao(1)    = qco * squeeze(magpcl);
phi2(1)      = squeeze(p0);

solver_jump_freq = NaN;             % first fminsearch failure (diagnostic)

% =========================================================================
%  Per-frequency harmonic-balance solve
% =========================================================================
for z = 1:N
    % Frequency response evaluation
    [mGc,  pGc ] = bode(Gc,  w(z));
    [mGac, pGac] = bode(Gac, w(z));
    magGc  = squeeze(mGc);   phGc  = squeeze(pGc);
    magGac = squeeze(mGac);  phGac = squeeze(pGac);

    % Harmonic-balance solve (warm-started from previous frequency)
    xo = [deltao(z), phi2(z)];
    [a, fval(z)] = fminsearch(@eqs, xo, optimset('Display','off'));

    % Describing-function evaluation at the converged operating point
    deltait = a(1);
    phi2it  = a(2);
    [magNit, phNit] = dfunction(w(z), R, deltait);

    A   = deltait * magGac * magNit / (K * qco);
    phi = phi2it + phGac + phNit*180/pi;

    % Closed-loop -> open-loop describing function (Nichols overlay)
    NoMag(z) = 20*log10( A / sqrt(1 - 2*A*cos(dtr*phi) + A^2) );
    NoPh(z)  = dtr*phi - atan2(-A*sin(dtr*phi), 1 - A*cos(dtr*phi));

    % Branch-tracking: restart on a known seed if the solver stalls
    if fval(z) > 1e-4
        if isnan(solver_jump_freq)
            solver_jump_freq = w(z);
        end
        deltao(z+1) = 24;          % values tuned for this specific aircraft
        phi2(z+1)   = -237;
    else
        deltao(z+1) = a(1);
        phi2(z+1)   = a(2);
    end
end

%% --------------------------------------------------------------------------
%  Linear OLOP onset frequency:  w * |pcl(jw)| = R     (FIX 9)
%  --------------------------------------------------------------------------
%  In the OLOP framework, the onset frequency is a property of the linear
%  closed loop, independent of the pilot amplitude qco. The rate-limit
%  line in the OLOP magnitude plot is  R/jw,  and it intersects |pcl(jw)|
%  at the onset frequency:
%
%        | pcl(jw_onset) | = R / w_onset    <=>    w * |pcl| = R
%
%  Linear interpolation gives a value finer than the frequency grid.
% --------------------------------------------------------------------------
[mag_pcl_v, ~] = bode(pcl, w);
mag_pcl_v      = squeeze(mag_pcl_v).';
rate_curve     = w .* mag_pcl_v;                    % w * |pcl(jw)|

idx_cross = find(rate_curve >= R, 1, 'first');
if isempty(idx_cross) || idx_cross == 1
    onset_freq = NaN;
    warning('No crossing of w*|pcl| = R found within the sweep.');
else
    pre = idx_cross - 1;
    onset_freq = interp1(rate_curve([pre idx_cross]), ...
                         w([pre idx_cross]), R, 'linear');
end

fprintf('OLOP onset frequency  (w*|pcl| = R)      : %.4g rad/s\n', onset_freq);
fprintf('Solver-jump frequency (first fval>1e-4)  : %.4g rad/s\n', solver_jump_freq);


%% Figure 1: Nichols chart aligned on (-180 deg, 0 dB)
H_dol   = (10.^(NoMag/20)) .* exp(1i * NoPh);
sys_dol = frd(H_dol(:), w, 'FrequencyUnit', 'rad/s');

H_lin   = squeeze(freqresp(Gc*Gac, w));
sys_lin = frd(H_lin, w, 'FrequencyUnit', 'rad/s');

opts                      = nicholsoptions;
opts.PhaseMatching        = 'on';
opts.PhaseMatchingFreq    = 1;
opts.PhaseMatchingValue   = -180;
opts.PhaseWrapping        = 'on';
opts.PhaseWrappingBranch  = -360;
opts.Grid                 = 'on';
opts.Title.String         = 'Nichols chart: linear loop vs. rate-limited OLOP DF';
opts.XLabel.String        = 'Open-loop phase';
opts.YLabel.String        = 'Open-loop magnitude';

figure(1);
hold on
h = nicholsplot(sys_lin, sys_dol, opts);

xline(-180, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1, 'HandleVisibility','off');
yline(   0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1, 'HandleVisibility','off');
plot(-180, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5, ...
     'MarkerFaceColor', 'r', 'HandleVisibility','off');

legend('Linear G_c G_{ac}', 'OLOP describing function', 'Location','best');
xlim([-300 -50]); ylim([-20 20]); grid on; hold off


%% Figure 2: Convergence diagnostic
figure(2);
semilogx(w, fval); grid on
xlabel('\omega [rad/s]');
ylabel('Residual f_{val}');
title('fminsearch residual vs. frequency');


%% Figure 3: w*|pcl| vs rate-limit line  -  shows the onset crossing
figure(3); clf
semilogx(w, 20*log10(rate_curve), 'b', ...
         w, 20*log10(R*ones(size(w))), 'k--');
hold on
xline(onset_freq, ':r', 'LineWidth', 1.2);
text(onset_freq, 20*log10(R)+1.5, ...
     sprintf('  \\omega_{onset} = %.4g rad/s', onset_freq), ...
     'Color', 'r');
grid on
xlabel('\omega [rad/s]')
ylabel('20 log_{10}( \omega \cdot |pcl| )   [dB]')
title('OLOP onset frequency: w \cdot |pcl(jw)| = R')
legend('w \cdot |pcl(jw)|', 'Rate limit R', 'Location','best')


%% THE LOCAL FUNCTIONS

function f = eqs(x)                             % harmonic-balance residual
    global w R dtr magGc magGac phGc phGac K qco z
    deltao = x(1);
    phi2   = x(2);

    [magN, phN] = dfunction(w(z), R, deltao);   % magN: gain, phN: rad

    re = deltao/magGc * cos(dtr*(phi2 - phGc)) + ...
         deltao*magGac*magN * cos(dtr*(phi2 + phGac) + phN) - K*qco;

    im = deltao/magGc * sin(dtr*(phi2 - phGc)) + ...
         deltao*magGac*magN * sin(dtr*(phi2 + phGac) + phN);

    f = re^2 + im^2;
end

function [magN, phN] = dfunction(freq, rate, inamp)
% Describing function of a rate-limiting element (Duda piecewise form).
%   magN  - dimensionless gain  (1 below onset; 4/(pi*x) far above)
%   phN   - phase in RADIANS    (0 below onset; -acos(pi/(2x)) far above)
    x = freq * inamp / rate;
    if x < 1                  % Region I
        magN = 1;
        phN  = 0;
    elseif x < 1.862          % Region II
        magN = polyval([ 0.2908 -1.4396  1.9232  0.2230 ], x);
        phN  = polyval([ 0.5280 -2.6213  3.5056 -1.4171 ], x);
    else                      % Region III
        magN = 4 / (x*pi);
        phN  = -acos(pi/(2*x));
    end
end
