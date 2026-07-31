%% fig8_v4.m
%  Fig. 8 style Nichols chart - negative inverse describing function
%  technique applied to the X-15, drawn with nicholsplot().
%
%  Two things make nicholsplot() behave like the hand-drawn figure:
%    1) opts.PhaseWrapping = 'on'  (+ PhaseWrappingBranch = -180 where
%       available) instead of the default 0...720 deg unwrapped axis;
%    2) both loci are handed over as FRD objects, and G(jw) is split into
%       continuous phase branches so the chart does not draw a chord across
%       the +/-180 deg wrap.
%
%  Rate limiter DF: Duda's three-region form, a function of the single
%  normalised argument alpha = w/w_onset, w_onset = R/A. The locus is built
%  on alpha, the crossing is located geometrically, and the frequency axis
%  of the -1/N FRD is then rescaled with the self-consistent limit-cycle
%  amplitude so that the displayed frequencies of the two curves agree at
%  the intersection - i.e. a true  G(jw_c) = -1/N(jw_c,u0)  crossing.
%
%  Requires: Control System Toolbox.
% -------------------------------------------------------------------------
clear; clc; close all;

%% ---------------------------- User switches -----------------------------
R                 = 15;             % actuator rate limit [deg/s]
GAIN_SCALE        = 1.0;            % airframe gain multiplier (0.1 / 10 tests)
PHASE_LIM         = [-180 -60];     % chart window, phase [deg]
MAG_LIM           = [  -5  15];     % chart window, amplitude [dB]
SHOW_NICHOLS_GRID = false;          % true -> M/N circles (ngrid)
SHOW_FREQ_TICKS   = true;           % frequency markers along G(jw)
FREQ_TICKS        = [2 3 4 5 8 12]; % [rad/s]
U0_FALLBACK       = 5;              % u0 [deg] used if no crossing is found
TARGET_W          = 2.8;            % reference LC frequency for the gain audit

%% --------- Linear plant : control law x aircraft dynamics ----------------
num_claw   = [5.21, -273.7855, -1425.456, -700.224];
den_claw   = [1, 21.36, 545.6, 605.7, 0];
Gs_claw    = tf(num_claw, den_claw);

num_ldynac = [-10.524, -16.8384, -0.6209, 0];
den_ldynac = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac  = tf(num_ldynac, den_ldynac);

Gs_ac      = GAIN_SCALE * Gs_claw * Gs_ldynac;   % linear OLTF  q / q_c

%% ------------------------- Locus 1 : G(jw) -------------------------------
w    = logspace(log10(0.5), log10(60), 4000);
Gjw  = reshape(squeeze(freqresp(Gs_ac, w)), 1, []);
phG  = wrap180(rad2deg(angle(Gjw)));
magG = 20*log10(abs(Gjw));

%% ------------------- Locus 2 : -1/N(jw,u0)  (Duda DF) --------------------
alpha       = linspace(1, 12, 4000);        % alpha = w/w_onset
[magN, phN] = rateLimiterDF(alpha);
negInvN     = -1 ./ (magN .* exp(1i*phN));
phD         = wrap180(rad2deg(angle(negInvN)));
magD        = 20*log10(abs(negInvN));

%% ------------------ Limit cycle : locus intersection ---------------------
S = polyIntersectLocal(phG, magG, w, phD, magD, alpha, 90);

fprintf('\n=== Describing-function limit cycle prediction (K_af x %.3g) ===\n', GAIN_SCALE);
if isempty(S)
    fprintf('  No intersection: DF predicts no limit cycle in the swept range.\n');
    U0 = U0_FALLBACK;
else
    for m = 1:numel(S)
        w_c = S(m).p1;  a_c = S(m).p2;  A_c = a_c*R/w_c;
        fprintf(['  #%d  w_c = %6.3f rad/s | phase = %7.2f deg | gain = %5.2f dB\n' ...
                 '      alpha_c = %5.3f  ->  A_LC = %6.2f deg , w_onset = %5.2f rad/s\n'], ...
                 m, w_c, S(m).x, S(m).y, a_c, A_c, w_c/a_c);
    end
    U0 = S(1).p2 * R / S(1).p1;             % self-consistent u0 [deg]
end
fprintf('  -1/N locus drawn for u0 = %.2f deg  (w_onset = %.2f rad/s)\n', U0, R/U0);

%% -------------------- FRD models handed to nicholsplot -------------------
% -1/N on a true frequency axis: alpha = w/w_onset, w_onset = R/u0
w_df   = alpha * R / U0;
sys_df = frd(negInvN, w_df);

% G(jw) split into continuous phase branches (no chord across the wrap)
segs = branchSegments(phG, 90);
args = {};
for s = 1:numel(segs)
    idx = segs{s};
    in  = phG(idx)  >= PHASE_LIM(1) & phG(idx)  <= PHASE_LIM(2) & ...
          magG(idx) >= MAG_LIM(1)   & magG(idx) <= MAG_LIM(2);
    if any(in) && numel(idx) > 1
        args = [args, {frd(Gjw(idx), w(idx)), 'k-'}];   %#ok<AGROW>
    end
end
args = [args, {sys_df, 'k--'}];

%% ------------------------- nicholsplot options ---------------------------
opts = nicholsoptions('cstprefs');
opts.PhaseUnits   = 'deg';
opts.MagUnits     = 'dB';
opts.FreqUnits    = 'rad/s';
opts.PhaseWrapping = 'on';                 % <-- key: no 0...720 unwrapping
try
    opts.PhaseWrappingBranch = -180;       % R2021a and later
catch
    % older releases already wrap to [-180,180)
end
opts.Grid         = ternary(SHOW_NICHOLS_GRID, 'on', 'off');
opts.XLim         = {PHASE_LIM};
opts.YLim         = {MAG_LIM};
opts.XLimMode     = {'manual'};
opts.YLimMode     = {'manual'};
opts.XLabel.String   = 'phase, deg';
opts.YLabel.String   = 'amplitude, dB';
opts.XLabel.FontSize = 11;
opts.YLabel.FontSize = 11;
opts.Title.String    = ['Negative inverse describing function technique - ' ...
                        'X-15, Nichols chart'];
opts.Title.FontSize  = 11;
opts.Title.FontWeight= 'normal';

%% ------------------------------- Plot ------------------------------------
figure('Color','w','Name','Fig. 8 - Negative inverse DF, X-15', ...
       'Position',[100 100 720 560]);
h  = nicholsplot(args{:}, opts);
ax = gca;

% line weights (nicholsplot LineSpec does not carry width)
set(findobj(ax,'Type','line','LineStyle','-' ),'LineWidth',1.6);
set(findobj(ax,'Type','line','LineStyle','--'),'LineWidth',1.4);

% ---- overlays (marker, labels, frequency ticks) --------------------------
try
    hold(ax,'on');

    plot(ax, PHASE_LIM, [0 0], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);

    for m = 1:numel(S)
        plot(ax, S(m).x, S(m).y, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 7);
        text(ax, S(m).x+3, S(m).y, sprintf('  limit cycle (%.1f rad/s)', S(m).p1), ...
             'FontSize', 10, 'VerticalAlignment','middle');
    end

    if SHOW_FREQ_TICKS
        phU = unwrap(angle(Gjw))*180/pi;
        for f = FREQ_TICKS
            pf = wrap180(interp1(w, phU, f));
            mf = interp1(w, magG, f);
            if pf >= PHASE_LIM(1) && pf <= PHASE_LIM(2) && ...
               mf >= MAG_LIM(1)   && mf <= MAG_LIM(2)
                plot(ax, pf, mf, 'k.', 'MarkerSize', 9);
                text(ax, pf+1.5, mf-0.6, sprintf('%g', f), 'FontSize', 8, ...
                     'Color', [0.3 0.3 0.3]);
            end
        end
    end

    text(ax, -120, 11.5, '-1/N(j\omega, u_0)', 'FontSize', 11);
    text(ax, -158,  1.5, 'G(j\omega)',         'FontSize', 11);

    hold(ax,'off');
catch ME
    warning('fig8:overlay', ...
        ['Annotation overlay skipped (%s). The chart itself is complete; ' ...
         'the crossing is listed in the command window.'], ME.message);
end

%% ------------------- Gain audit for the factor-of-10 check ---------------
phT  = wrap180(interp1(w, unwrap(angle(Gjw))*180/pi, TARGET_W));
magT = interp1(w, magG, TARGET_W);
[phD_s, is] = sort(phD);  magD_s = magD(is);
[phD_s, iu] = unique(phD_s);  magD_s = magD_s(iu);
magD_T = interp1(phD_s, magD_s, phT, 'linear', NaN);
if ~isnan(magD_T)
    k_req = 10^((magD_T - magT)/20);
    fprintf(['\n  Gain audit: to place the limit cycle at %.2f rad/s the linear\n' ...
             '  loop gain must be scaled by %.3f (%.2f dB) w.r.t. the current model.\n'], ...
             TARGET_W, k_req, 20*log10(k_req));
end

%% ============================ local functions ============================
function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function p = wrap180(p)
p = mod(p + 180, 360) - 180;                 % -> (-180, 180]
end

function [m, ph] = rateLimiterDF(alpha)
% Duda's three-region rate-limiter describing function.
% alpha = w/w_onset ; m [-] linear magnitude ; ph [rad] phase.
m  = zeros(size(alpha));
ph = zeros(size(alpha));
for k = 1:numel(alpha)
    a = alpha(k);
    if a < 1                                  % Region I  : no saturation
        m(k)  = 1;              ph(k) = 0;
    elseif a < 1.862                          % Region II : transition
        m(k)  = 0.2908*a^3 - 1.4396*a^2 + 1.9232*a + 0.223;
        ph(k) = 0.5280*a^3 - 2.6213*a^2 + 3.5056*a - 1.4171;
    else                                      % Region III: full saturation
        wbar  = 1/a;
        m(k)  = 4*wbar/pi;      ph(k) = -acos(min(1, pi*wbar/2));
    end
end
end

function segs = branchSegments(ph, thr)
% Split a wrapped phase vector into index ranges free of +/-180 deg jumps.
k     = [0, find(abs(diff(ph)) > thr), numel(ph)];
segs  = cell(1, numel(k)-1);
for m = 1:numel(k)-1
    segs{m} = (k(m)+1):k(m+1);
end
end

function S = polyIntersectLocal(x1, y1, p1, x2, y2, p2, jumpTol)
% Intersections of two polylines, with the interpolated curve parameters
% (frequency on curve 1, alpha on curve 2) returned at each crossing.
S = struct('x',{},'y',{},'p1',{},'p2',{});
dx2 = diff(x2);  dy2 = diff(y2);
ok2 = abs(dx2) < jumpTol;
for i = 1:numel(x1)-1
    rx = x1(i+1) - x1(i);  ry = y1(i+1) - y1(i);
    if abs(rx) >= jumpTol, continue; end
    d = rx*dy2 - ry*dx2;
    d(abs(d) < eps) = NaN;
    qx = x2(1:end-1) - x1(i);
    qy = y2(1:end-1) - y1(i);
    t  = (qx.*dy2 - qy.*dx2) ./ d;
    u  = (qx.*ry  - qy.*rx ) ./ d;
    j  = find(ok2 & t >= 0 & t <= 1 & u >= 0 & u <= 1);
    for m = 1:numel(j)
        k = j(m);
        S(end+1).x = x1(i) + t(k)*rx;                     %#ok<AGROW>
        S(end).y   = y1(i) + t(k)*ry;
        S(end).p1  = p1(i) + t(k)*(p1(i+1) - p1(i));
        S(end).p2  = p2(k) + u(k)*(p2(k+1) - p2(k));
    end
end
if numel(S) > 1                                % merge duplicates
    keep = true(1, numel(S));
    for a = 2:numel(S)
        for b = 1:a-1
            if keep(b) && hypot(S(a).x - S(b).x, S(a).y - S(b).y) < 1e-3
                keep(a) = false; break
            end
        end
    end
    S = S(keep);
end
end