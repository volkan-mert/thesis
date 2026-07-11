%% ------------------------------------------------------------------------
%  Reproduction of Figure 2-8  "Time Domain Analysis"
%  Gilbreath (2001), "Prediction of Pilot-Induced Oscillations (PIO) Due to
%  Actuator Rate Limiting Using the Open-Loop Onset Point (OLOP) Criterion"
%
%  Nonlinear time-domain simulation of the closed-loop system of Figure 2-3
%  (highly augmented aircraft with a rate limiter inside the loop). The model
%  data are taken verbatim from the AIAA-95-3204-CP example coded in the
%  thesis Appendix B M-files:
%
%       qc --->(K)--->( + )--->[ Gc ]--->[ Rate Limiter (R) ]--->[ Gac ]---> q
%                       ^-                                                 |
%                       |_________________________________________________|
%
%  The Nichols/OLOP analysis predicts an onset frequency ~5.1 rad/sec.
%  This script confirms it in the time domain:
%     omega = 5.0 rad/sec  -> stable, bounded response
%     omega = 5.3 rad/sec  -> rate limiter activates, loop goes unstable
%                             after ~3 s (as described in the thesis text).
%
%  Requires: Control System Toolbox (tf2ss).  No Simulink needed.
%  Author of reproduction: (drop-in script)
%% ------------------------------------------------------------------------
clear; clc; close all;

%% ---------------------- System definition (Appendix B) ------------------
K   = 13.68;      % command / pilot-station gain
R   = 60;         % actuator rate limit                       [deg/sec]
qco = 1.1;        % pitch-rate command amplitude              [deg/sec]

% Control law  Gc = 5.21 (s-57.36)(s+4.26)(s+0.55)
%                   -------------------------------------------------------
%                   (s^2 + 2*0.442*22.85 s + 22.85^2) * s * (s+1.16)
numGc = 5.21 * conv(conv([1 -57.36],[1 4.26]),[1 0.55]);
denGc = conv(conv([1 2*0.442*22.85 22.85^2],[1 0]),[1 1.16]);
tf_Gc = tf(numGc,denGc)
% Aircraft dynamics  Gac = -10.524 (s+1.562)(s+0.038) s
%                          --------------------------------------------------
%                          (s^2 + 2*0.212*0.088 s + 0.088^2)(s+3.75)(s-1.44)
numGac = -10.524 * conv(conv([1 1.562],[1 0.038]),[1 0]);
denGac = conv(conv([1 2*0.212*0.088 0.088^2],[1 3.75]),[1 -1.44]);
tfGac = tf(numGac,denGac)
% State-space realizations (controllable canonical form)
[Ac,Bc,Cc,Dc] = tf2ss(numGc ,denGc );   nc = size(Ac,1);
[Aa,Ba,Ca,Da] = tf2ss(numGac,denGac);    na = size(Aa,1);

%% ---------------------- Simulation settings -----------------------------
Kh    = 2000;          % high gain used to realise an *ideal* rate limiter
                       % via a saturated integrator  d(delta)/dt = sat(Kh(dc-delta),R)
tmax  = 5;             % [sec]
tspan = linspace(0,tmax,3000);
x0    = zeros(nc+na+1,1);       % [ Gc states ; Gac states ; delta ]

omegas = [5.0 5.3];             % rad/sec  (stable / unstable case)
opts   = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',2e-3);

%% ---------------------- Run + plot (2x2, as in Fig 2-8) -----------------
figure('Color','w','Position',[80 80 1050 680]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:numel(omegas)
    w = omegas(j);

    % integrate the nonlinear closed loop
    [t,X] = ode15s(@(t,x) loopODE(t,x,w,nc,na,K,qco,R,Kh,Ac,Bc,Cc,Dc,Aa,Ba,Ca,Da), tspan, x0, opts);

    % reconstruct signals
    xc = X(:,1:nc).';  xa = X(:,nc+1:nc+na).';  delta = X(:,end).';
    q   = (Ca*xa + Da*delta).';                  % pitch rate           [deg/s]
    qc  = qco*sin(w*t);                          % pilot command        [deg/s]
    e   = K*qc - q;                              % loop error
    dc  = (Cc*xc).' + Dc*e;                      % rate-limiter input   [deg]
    dlt = delta.';                               % rate-limited output  [deg]

    % ---- top row: pitch rate ----
    nexttile(j);
    plot(t,q,'b-','LineWidth',1.4); hold on;
    plot(t,K*qc,'r--','LineWidth',1.0);
    grid on; box on;
    title(sprintf('\\omega = %.1f rad/sec',w),'FontWeight','bold');
    ylabel('pitch rate  [deg/s]');
    legend({'q  (output)','K\cdotq_c  (command)'},'Location','best','FontSize',8);
    
    if w > 5.1, ylim([-400 400]); end          % unstable case: show divergence

    % ---- bottom row: rate limiter element ----
    nexttile(j+2);
    plot(t,dc ,'m-','LineWidth',1.0); hold on;
    plot(t,dlt,'k-','LineWidth',1.4);
    grid on; box on;
    xlabel('Time (sec)'); ylabel('RLE  [deg]');
    legend({'\delta_c  (commanded)','\delta  (rate limited)'}, ...
           'Location','best','FontSize',8);
    
    if w > 5.1, ylim([-400 400]); end

    % annotate where the limiter saturates (|d(delta)/dt| pinned at R)
    active = abs(Kh*(dc-dlt)) >= R;

    if any(active)
        yl = ylim;
        text(t(find(active,1,'last')), 0.75*yl(2), ' rate limiting', 'Color',[0.3 0.3 0.3],'FontAngle','italic','FontSize',9);
    end
end

title(tl,'Figure 2-8 (reproduced):  Time Domain Analysis', 'FontWeight','bold','FontSize',12);

%% ---------------------- Local ODE function ------------------------------
function dx = loopODE(t,x,w,nc,na,K,qco,R,Kh,Ac,Bc,Cc,Dc,Aa,Ba,Ca,Da)
    xc    = x(1:nc);
    xa    = x(nc+1:nc+na);
    delta = x(end);

    q   = Ca*xa + Da*delta;          % plant output (pitch rate)
    qc  = qco*sin(w*t);              % sinusoidal pilot command
    e   = K*qc - q;                  % summing junction
    dc  = Cc*xc + Dc*e;              % control-law output = rate-limiter input

    % ideal rate limiter as a saturated high-gain integrator
    ddelta = max(-R, min(R, Kh*(dc - delta)));

    dxc = Ac*xc + Bc*e;
    dxa = Aa*xa + Ba*delta;
    dx  = [dxc; dxa; ddelta];
end
