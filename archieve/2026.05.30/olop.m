clear; clc; close all
format short

% Forward-loop gain
K = 13.68;
% Controller (pilot error -> RLE input)
Gc = tf( 5.21 * conv([1 -57.36], conv([1 4.26], [1 0.55])), conv([1 2*0.442*22.85 22.85^2], conv([1 0], [1 1.16])) );
% Actuator dynamics (RLE output -> theta)
Gac = tf(-10.524 * conv([1 1.562], conv([1 0.038], [1 0])), conv([1 2*0.212*0.088 0.088^2], conv([1 3.75], [1 -1.44])) );

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
Yc = ss( K*Gc*Gac / (1 + Gc*Gac) );

% ----- Pure-gain pilot search ----------------------------------------------
warning('off', 'Control:analysis:MarginUnstable')
Kg = (0.1*2.5):0.01:10;
for i = 1:length(Kg)
    [gm, pm, wcg, wcp] = margin(Kg(i)*Yc);
    if tag == 4, break, end
    if -180 + pm > xover
        Pilot_gain  = Kg(i);
        break
    end
end
warning('on', 'Control:analysis:MarginUnstable')

% ----- Neal-Smith modified pilot model -------------------------------------
Yp  = Kp * tf([5 1],[1 0]) * tf([Tp1 1],[Tp2 1]);
set(Yp, 'InputDelay', 0.25);
Ypp = pade(Yp, 2);                      % Pade approximation for ss/series ops

% ----- Closed-loop RLE input / pilot_stick ---------------------------------
P_RLE = ss( K*Gc / (1 + Gc*Gac) );
[magP, phP] = bode(P_RLE, w);

% ===========================================================================
%                    Loop over rate limits
% ===========================================================================
figure(3)
hold on

for n = 1:length(R)
    % Rate-limit line for the frequency plot:  R / (j*omega)
    [magR, phR] = bode(tf(R(n), [1 0]), w);
    
    % Closed-loop onset frequency: |Po * P_RLE(jw)| meets the rate-limit line
    iter = 100;
    for i = 1:length(w)
        dif(i) = abs(magR(i) - Po*magP(i));
        if dif(i) <= iter
            iter = dif(i);
            k    = i;
        end
    end
    
    % ----- OLOP transfer function -----------------------------------------
    if tag == 4
        Pilot = Ypp;                    % Neal-Smith with Pade-approx delay
    else
        Pilot = tf(Pilot_gain, 1);      % pure-gain pilot
    end
    OLOP = Pilot * Yc;                  % pilot * controlled element
    
    [magN, phN] = nichols(OLOP, w);
    magN = squeeze(magN).';
    phN  = squeeze(phN ).';
    
    % ----- OLOP point + locus on the Nichols chart ------------------------
    r  = 150;                           % plotted range around onset index
    
    markerStyles = {'o','d','o','p','s','d','s'};
    markerFaces  = {'k','k','none','k','k','none','none'};
    idx = min(n, numel(markerStyles));
    
    plot(phN(k), 20*log10(magN(k)), markerStyles{idx}, 'MarkerFaceColor', markerFaces{idx}, 'MarkerSize', 7);
    
    iLo = max(1,        k - r);
    iHi = min(numel(w), k + r);
    plot(phN(iLo:iHi), 20*log10(magN(iLo:iHi)), [-360 0], [0 0], '-k', [-180 -180], [-50 50], '-k')
end

% ----- Stability Boundary and Formatting ----------------------------------
v1 = [-60 -90 -100 -120 -140 -160 -180];
v2 = [13.5 7.5 5.5 2.5 1.1 0 0];    % OLOP stability boundary

axis([-200, -60, -10, 20]);
plot(v1, v2, '-k', 'LineWidth', 1.5)    
grid on

xlabel('Open-loop phase [deg]')
ylabel('Open-loop magnitude [dB]')
title('OLOP Nichols chart with stability boundary')