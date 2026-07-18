%***************************************** olop.m ***************************
%  Open-loop Onset Point (OLOP) analysis for the aircraft system with a
%  rate-limiting element.
%
%  Aircraft model: AIAA-95-3204-CP (Duda / Hanke / Gilbreath).
%
%=============================================================================
%  CORRECTIONS APPLIED TO THE ORIGINAL gilbreathOLOP.m by CLAUDE AI
%=============================================================================
%
%  FIX 1 - LINMOD removed (the original Simulink models nt33v2Dlin.mdl and
%          nt33v2dolop.mdl are not available). The two linmod calls are
%          replaced with analytical closed-loop transfer functions:
%
%            Yc    = ss(A ,B ,C (4,:),D (4,:))  ->  K*Gc*Gac / (1 + Gc*Gac)
%            P_RLE = ss(A ,B ,C (6,:),D (6,:))  ->  K*Gc      / (1 + Gc*Gac)
%            OLOP  = ss(AA,BB,-CC    ,-DD)      ->  Pilot * Yc
%
%          The "-CC,-DD" sign flip in the original was a linmod
%          feedback-sign convention adjustment; the analytical form needs
%          no flip because Yc already carries the correct sign.
%
%  FIX 2 - lin_abcd.m inlined. The aircraft model constants K, Gc, Gac
%          are defined directly below this header instead of being loaded
%          from an external script. The script is now self-contained.
%
%  FIX 3 - ngridneal embedded as a local function at the bottom of this
%          file. The original referenced an ngridneal.m that ships with
%          the Gilbreath OLOP package but is not in standard MATLAB.
%          The local copy draws the standard Nichols grid plus the +3 dB
%          M-contour (Neal-Smith Level-1/2 closed-loop resonance bound).
%
%  FIX 4 - Margin-warning suppression bracket around the Kg pilot-gain
%          sweep. margin(Kg(i)*Yc) emits "closed-loop unstable" warnings
%          at high gains; the bracket silences only that specific warning
%          and restores it after the sweep.
%
%  FIX 5 - nichols(OLOP,w) output squeeze. Like bode, nichols returns
%          1xNxNw 3-D arrays for SISO inputs. magN and phN are squeezed
%          to row vectors so phN(k), phN(iLo:iHi) etc. work in 1-D.
%
%  FIX 6 - Locus-slice index clipping. The original phN(1,k-r:r+k) can
%          go out of bounds when k < r+1 or k+r > numel(w). The indices
%          are now clamped:
%              iLo = max(1,        k - r)
%              iHi = min(numel(w), k + r)
%
%  FIX 7 - Marker-style ladder replaced with cell arrays indexed by n
%          (markerStyles, markerFaces). Identical visual result, far
%          shorter source, and trivial to extend past 7 rate-limit cases.
%
%  FIX 8 - Cosmetic: xlabel / ylabel / title on figure(4); consistent
%          section banners and indentation. No numerical impact.
%
%          visdiff('gilbreathOLOP.m','gilbreathOLOP_original.m')
%
%=============================================================================
%  Analytical relations used in FIX 1:
%
%       Yc    = theta / pilot_stick   = K * Gc * Gac / (1 + Gc * Gac)
%       P_RLE = RLE_in / pilot_stick  = K * Gc      / (1 + Gc * Gac)
%       OLOP  = Pilot * Yc
%=============================================================================

clear; clc; close all
format short

% ===========================================================================
%  FIX 2: Aircraft model (AIAA-95-3204-CP)  -  inlined from lin_abcd.m
% ===========================================================================
%  Loop topology (rate limiter shown for context, identity in linear model)
%
%       qco --> [ K ] --(+)--> [ Gc ] --> delta_c --[ RLE ]--> delta
%                       ^                                         |
%                       |                                         v
%                       '----------- theta <-------[ Gac ]<--------'

% Forward-loop gain
K = 13.68;

% Controller (pilot error -> RLE input)
Gc = tf( 5.21 * conv([1 -57.36], conv([1 4.26], [1 0.55])), conv([1 2*0.442*22.85 22.85^2], conv([1 0], [1 1.16])) );

% Actuator dynamics (RLE output -> theta)
Gac = tf(-10.524 * conv([1 1.562], conv([1 0.038], [1 0])), conv([1 2*0.212*0.088 0.088^2], conv([1 3.75], [1 -1.44])) );

% ===========================================================================
%  OLOP analysis
% ===========================================================================

w = logspace(-2, 2, 1000);

% ----- User inputs ----------------------------------------------------------
R = input(' Enter Vector of Rate Limits to evaluate [default = 60]: ');
if isempty(R), R = 60; end

Per = input(' Enter Stick Input Amplitude in percentage (range 0-100) [default = 100%]: ');
if isempty(Per), Per = 100; end

Po = Per/100 * 3.6;
if Po < 0.8                            % non-linear stick gradient
    Po = 2.5 * Po;
else
    Po = 2.5 * (2.86*Po - 1.48);
end

tag = input(' Choose desired pilot model (1=Low 2=Med 3=High 4=Neal-Smith) [default Neal-Smith]: ');
if isempty(tag)
    tag = 4;
elseif tag == 1, xover = -90;
elseif tag == 2, xover = -110;
elseif tag == 3, xover = -130;
end

% ----- Pilot model parameters (Neal-Smith) ---------------------------------
Kp  = -0.05 * 2.5;
Tp1 = 0.06;
Tp2 = 0.01;

% ----- Controlled element (theta / pilot_stick, closed loop) ---------------
% FIX 1: was  [A,B,C,D] = LINMOD('nt33v2Dlin'); Yc = ss(A,B,C(4,:),D(4,:));
Yc = ss( K*Gc*Gac / (1 + Gc*Gac) );

% ----- Pure-gain pilot search ----------------------------------------------
% FIX 4: silence the "closed-loop unstable" margin warnings during the sweep
warning('off', 'Control:analysis:MarginUnstable')

Kg = (0.1*2.5):0.01:10;
for i = 1:length(Kg)
    [gm, pm, wcg, wcp] = margin(Kg(i)*Yc);
    if tag == 4, break, end
    if -180 + pm > xover
        xover_angle = -180 + pm
        Pilot_gain  = Kg(i)
        break
    end
end

warning('on', 'Control:analysis:MarginUnstable')

figure(1), ngridneal, nichols(Kg(i)*Yc, w)        % FIX 3: ngridneal is local
grid on
% ----- Neal-Smith modified pilot model -------------------------------------
bw  = 3.5;                              % bandwidth requirement
Yp  = Kp * tf([5 1],[1 0]) * tf([Tp1 1],[Tp2 1]);
set(Yp, 'InputDelay', 0.25);
Ypp = pade(Yp, 2);                      % Pade approximation for ss/series ops

figure(2), ngridneal, nichols(Yp*Yc, {0.1, 3*bw}) % FIX 3: ngridneal is local
grid on
% ----- Closed-loop RLE input / pilot_stick ---------------------------------
% FIX 1: was  P_RLE = ss(A,B,C(6,:),D(6,:));
P_RLE = ss( K*Gc / (1 + Gc*Gac) );
[magP, phP] = bode(P_RLE, w);

% ===========================================================================
%                    Loop over rate limits
% ===========================================================================
for n = 1:length(R)

    % Rate-limit line for the frequency plot:  R / (j*omega)
    [magR, phR] = bode(tf(R(n), [1 0]), w);

    figure(3)
    semilogx(w, 20*log10(Po*magP(1,:)), w, 20*log10(magR(1,:)), 'k:')
    hold on, axis([0.01, 10, 10, 60])
    title('Determine Closed-loop Onset Frequency')

    % Closed-loop onset frequency: |Po * P_RLE(jw)| meets the rate-limit line
    iter = 100;
    for i = 1:length(w)
        dif(i) = abs(magR(i) - Po*magP(i));
        if dif(i) <= iter
            iter = dif(i);
            k    = i;
        end
    end
    w_onset = w(k)
    Rate = R(n)

    % ----- OLOP transfer function -----------------------------------------
    % FIX 1: was  [AA,BB,CC,DD] = linmod('nt33v2dolop');
    %             OLOP = ss(AA,BB,-CC,-DD);
    if tag == 4
        Pilot = Ypp;                    % Neal-Smith with Pade-approx delay
    else
        Pilot = tf(Pilot_gain, 1);      % pure-gain pilot
    end
    OLOP = Pilot * Yc;                  % pilot * controlled element

    % FIX 5: squeeze 3-D nichols output to row vectors so 1-D indexing works
    [magN, phN] = nichols(OLOP, w);
    magN = squeeze(magN).';
    phN  = squeeze(phN ).';

    % ----- OLOP point + locus on the Nichols chart ------------------------
    r  = 150;                           % plotted range around onset index
    v1 = [-60 -90 -100 -120 -140 -160 -180];
    v2 = [13.5 7.5 5.5 2.5 1.1 0 0];    % OLOP stability boundary

    % FIX 7: marker ladder collapsed into indexed cell arrays
    markerStyles = {'o','d','o','p','s','d','s'};
    markerFaces  = {'k','k','none','k','k','none','none'};
    idx = min(n, numel(markerStyles));

    figure(4)
    plot(phN(k), 20*log10(magN(k)), markerStyles{idx}, 'MarkerFaceColor', markerFaces{idx}, 'MarkerSize', 7);
    grid on
    hold on

    % FIX 6: clamp the locus slice so it never indexes out of bounds
    iLo = max(1,        k - r);
    iHi = min(numel(w), k + r);
    plot(phN(iLo:iHi), 20*log10(magN(iLo:iHi)), [-360 0], [0 0], '-k', [-180 -180], [-50 50], '-k')

    axis([-200, -60, -10, 20]);
    plot(v1, v2, '-k', 'LineWidth', 1.5)    % OLOP stability boundary
    grid on
end

% FIX 8: figure-4 labels
xlabel('Open-loop phase [deg]')
ylabel('Open-loop magnitude [dB]')
title('OLOP Nichols chart with stability boundary')


%% =========================================================================
%  FIX 3: ngridneal embedded as a local function (was a missing helper)
%  =========================================================================

function ngridneal
%NGRIDNEAL  Nichols grid for Neal-Smith / OLOP pilot-vehicle analysis.
%
%   Stand-in for the (missing) ngridneal helper from the original
%   Gilbreath OLOP package. Draws:
%     (1) the standard Nichols grid (closed-loop M and N contours)
%     (2) a +3 dB M-contour highlighted in red - the Neal-Smith
%         Level-1 / Level-2 closed-loop resonant-peak boundary used in
%         pilot handling-qualities analysis.

    ngrid                                             

    hold on
    ax = gca;

    % Highlight the +3 dB M-circle.
    % Parametric trace of |T| = M in the L = T/(1-T) plane:
    M_dB = 3;
    M    = 10^(M_dB/20);
    phi  = linspace(eps, 2*pi - eps, 720);
    T    = M * exp(1i*phi);
    L    = T ./ (1 - T);

    magL_dB = 20*log10(abs(L));
    phaseL  = unwrap(angle(L)) * 180/pi;
    phaseL  = mod(phaseL, 360) - 360;   % wrap to [-360, 0]

    plot(ax, phaseL, magL_dB, 'r-', 'LineWidth', 1.4, 'HandleVisibility','off');
    text(ax, -180, M_dB + 1, sprintf('  M = %+g dB (NS bound)', M_dB), 'Color', 'r', 'FontSize', 8, 'HandleVisibility','off');
    grid on

    hold off
end